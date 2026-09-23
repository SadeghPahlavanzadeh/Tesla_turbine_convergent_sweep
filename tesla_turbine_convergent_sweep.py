import numpy as np
import pandas as pd
import ansys.fluent.core as pyfluent

# ==============================================================================
# 1. LAUNCH FLUENT SESSION
# ==============================================================================
session = pyfluent.launch_fluent(
    precision="double",
    processor_count=4,
    mode="solver",
    show_gui=False
)
print("Fluent session launched successfully.")

# ==============================================================================
# 2. READ BASELINE CASE & SETUP TURBULENCE MODEL
# ==============================================================================
CASE_FILE = "tesla_turbine.cas.h5"
print(f"Loading baseline case: {CASE_FILE}")
session.tui.file.read_case(CASE_FILE)

# Set k-omega SST turbulence model
session.setup.models.viscous.model = "k-omega"
session.setup.models.viscous.k_omega_model = "sst"
print("Baseline case loaded and k-omega SST model active.")

# ==============================================================================
# 3. HELPER FUNCTIONS & ADVANCED CONVERGENCE CRITERIA
# ==============================================================================
P_ATM = 101325.0  # Atmospheric pressure in Pascals

def update_boundary_conditions(session, p_inlet_abs_bar, rpm_disk):
    """
    Updates the plenum inlet absolute pressure (converted to gauge for Fluent)
    and disk wall angular velocity in rad/s.
    """
    p_inlet_gauge_pa = (p_inlet_abs_bar * 1e5) - P_ATM
    
    # Set Pressure Inlet Boundary
    session.setup.boundary_conditions.pressure_inlet["plenum_inlet"].gauge_total_pressure.value = p_inlet_gauge_pa
    
    # Set Disk Rotational Speed (rad/s)
    omega_rad_s = rpm_disk * (2.0 * np.pi / 60.0)
    session.setup.boundary_conditions.wall["disk_wall"].angular_velocity.value = omega_rad_s
    
    return omega_rad_s

def check_and_refine_mesh(session):
    """
    Evaluates y+ and minimum domain pressure from the completed solution,
    executing boundary layer refinement if thresholds are exceeded.
    """
    max_yplus = session.solution.reports.surface.maximum(
        surface_names=["disk_wall"], expression="y-plus"
    )
    min_gauge_p = session.solution.reports.volume.minimum(
        zone_names=["fluid_domain"], expression="pressure"
    )
    
    refined = False
    
    # Rule 1: High expansion / low gauge pressure condition
    if min_gauge_p < -3000.0:
        print(f"      [Refinement Triggered] Min Gauge Pressure ({min_gauge_p:.1f} Pa) < -3000 Pa. Refining first 4 prism layers...")
        session.tui.adapt.boundary_layer("disk_wall", "yes", "4")
        refined = True
        
    # Rule 2: Boundary layer resolution condition
    elif max_yplus > 1.0:
        print(f"      [Refinement Triggered] Max Y+ ({max_yplus:.3f}) > 1.0. Refining first prism layer...")
        session.tui.adapt.boundary_layer("disk_wall", "yes", "1")
        refined = True

    return refined, max_yplus, min_gauge_p

def check_mass_balance(session, inlet_zone="plenum_inlet", outlet_zone="outlet", threshold_pct=0.1):
    """
    1. Checks mass balance conservation: relative mass flow imbalance between inlet and outlet.
    """
    try:
        m_in = abs(session.solution.reports.surface.mass_flow_rate(surface_names=[inlet_zone]))
        net_imbalance = abs(session.solution.reports.surface.mass_flow_rate(surface_names=[inlet_zone, outlet_zone]))
        if m_in > 0:
            imbalance_pct = (net_imbalance / m_in) * 100.0
            return imbalance_pct < threshold_pct, imbalance_pct
    except Exception:
        pass
    return False, 100.0

def monitor_torque_stability(session, wall_name="disk_wall", history_list=None, window_size=5, tolerance=0.001):
    """
    2. Checks torque generation stability as the integral of shear stress and pressure moments on disks.
    """
    if history_list is None:
        history_list = []
        
    current_torque = session.solution.reports.surface.moment(
        surface_names=[wall_name], moment_center=[0, 0, 0], axis=[0, 0, 1]
    )
    history_list.append(current_torque)
    
    if len(history_list) > window_size:
        history_list.pop(0)
        
    if len(history_list) == window_size:
        torque_range = max(history_list) - min(history_list)
        mean_torque = abs(np.mean(history_list))
        if mean_torque > 0:
            relative_variation = torque_range / mean_torque
            return relative_variation < tolerance, current_torque, history_list
            
    return False, current_torque, history_list

def track_min_pressure_fluctuation(session, iteration_count, min_p_tracker):
    """
    3. After 5000 iterations, captures local minimum pressure and monitors 
       its fluctuation until pressure change is lower than 100 Pa.
    """
    current_min_p = session.solution.reports.volume.minimum(
        zone_names=["fluid_domain"], expression="pressure"
    )
    
    if iteration_count >= 5000 and not min_p_tracker["active"]:
        min_p_tracker["active"] = True
        min_p_tracker["last_val"] = current_min_p
        print(f"      [Min P Tracker] Reached 5000 iterations. Initial Min Pressure: {current_min_p:.2f} Pa. Tracking delta (< 100 Pa)...")
        return False

    if min_p_tracker["active"]:
        delta_p = abs(current_min_p - min_p_tracker["last_val"])
        min_p_tracker["last_val"] = current_min_p
        print(f"      [Min P Tracker] Current Min P: {current_min_p:.2f} Pa | Delta: {delta_p:.2f} Pa")
        if delta_p < 100.0:
            return True
            
    return False

def solve_with_advanced_convergence(session, max_total_iterations=10000, chunk_size=100, max_passes=2):
    """
    Executes iterative solver chunks coupled with mesh adaptation and advanced 
    physical/numerical convergence checks.
    """
    session.solution.initialization.hybrid_initialize()
    
    pass_cnt = 0
    refined, final_yplus, min_gauge_p = check_and_refine_mesh(session)
    
    # Outer adaptive loop for mesh quality
    while pass_cnt < max_passes:
        if refined:
            session.solution.initialization.hybrid_initialize()
            
        total_iters = 0
        torque_history = []
        min_p_tracker = {"active": False, "last_val": 0.0}
        converged = False
        
        # Inner chunk-based convergence loop
        while total_iters < max_total_iterations:
            session.solution.run_calculation.iterative(number_of_iterations=chunk_size)
            total_iters += chunk_size
            
            # Evaluate Convergence Criteria
            mass_ok, imbalance = check_mass_balance(session)
            torque_ok, current_torque, torque_history = monitor_torque_stability(session, history_list=torque_history)
            min_p_ok = track_min_pressure_fluctuation(session, total_iters, min_p_tracker)
            
            # Print intermediate status
            print(f"   [Iter {total_iters}] Mass Imbalance: {imbalance:.3f}% | Torque: {current_torque:.5f} Nm | Min P Active: {min_p_tracker['active']}")
            
            # Exit criteria: Must exceed 5000 iterations for min pressure stability check, plus stable mass balance and torque
            if total_iters >= 5000 and mass_ok and torque_ok and min_p_ok:
                print(f"   --> Advanced convergence criteria satisfied at iteration {total_iters}.")
                converged = True
                break
        
        # Post-solve mesh check
        refined, final_yplus, min_gauge_p = check_and_refine_mesh(session)
        if not refined or pass_cnt >= max_passes - 1:
            break
        pass_cnt += 1

    return final_yplus, min_gauge_p

def calculate_performance_metrics(session, omega_rad_s, p_inlet_abs_bar):
    """
    Calculates Euler generated power (P = Torque * omega) and total-to-static
    isentropic efficiency based on total inlet pressure.
    """
    torque = session.solution.reports.surface.moment(
        surface_names=["disk_wall"], moment_center=[0, 0, 0], axis=[0, 0, 1]
    )
    power_euler = torque * omega_rad_s
    m_dot = abs(session.solution.reports.surface.mass_flow_rate(surface_names=["plenum_inlet"]))
    
    gamma = 1.4
    cp = 1004.5   # J/(kg·K)
    t01 = 298.15  # K
    p01_pa = p_inlet_abs_bar * 1e5
    p2_pa = P_ATM
    
    p_isentropic = m_dot * cp * t01 * (1.0 - (p2_pa / p01_pa) ** ((gamma - 1.0) / gamma))
    
    if p_isentropic > 0:
        efficiency_pct = (power_euler / p_isentropic) * 100.0
    else:
        efficiency_pct = 0.0
        
    return power_euler, torque, m_dot, efficiency_pct

# ==============================================================================
# 4. EXECUTE PARAMETRIC SWEEP
# ==============================================================================
p_inlet_levels_bar = np.round(np.arange(2.0, 3.0 + 0.05, 0.1), 1)
rpm_levels = np.arange(12000, 25000 + 1, 500)

results_table = []
run_count = 1

print(f"\nStarting 297-point parametric sweep with advanced convergence checks...")

for p_bar in p_inlet_levels_bar:
    print(f"\n==================================================")
    print(f"--- INLET PRESSURE: {p_bar:.1f} bar (abs) ---")
    print(f"==================================================")
    
    for rpm in rpm_levels:
        run_id = f"DP_{run_count:03d}_P{p_bar:.1f}bar_{rpm}RPM"
        print(f"\nRunning {run_id}: Inlet P = {p_bar:.1f} bar | Speed = {rpm} RPM")
        
        # 1. Update boundary conditions
        omega = update_boundary_conditions(session, p_bar, rpm)
        
        # 2. Run solver with advanced convergence loops & mesh refinement
        final_yplus, min_gauge_p = solve_with_advanced_convergence(session)
        
        # 3. Calculate Performance Metrics
        power_w, torque_nm, m_dot_kg_s, eta = calculate_performance_metrics(session, omega, p_bar)
        
        # 4. Store results
        results_table.append({
            "Run_ID": run_id,
            "Inlet_P_Abs_bar": p_bar,
            "Rotational_Speed_RPM": rpm,
            "Euler_Power_W": power_w,
            "Torque_Nm": torque_nm,
            "Mass_Flow_Rate_kg_s": m_dot_kg_s,
            "Efficiency_Pct": eta,
            "Final_Max_YPlus": final_yplus,
            "Min_Gauge_Pressure_Pa": min_gauge_p
        })
        
        print(f"   -> Power: {power_w:.2f} W | Efficiency: {eta:.2f}% | Y+: {final_yplus:.3f}")
        run_count += 1

# Exit session upon completion
session.exit()
print("\nParametric sweep completed successfully!")

# ==============================================================================
# 5. EXPORT AND DISPLAY RESULTS
# ==============================================================================
df = pd.DataFrame(results_table)
df.to_csv("tesla_turbine_advanced_parametric_results.csv", index=False)
print("Results successfully exported to tesla_turbine_advanced_parametric_results.csv")