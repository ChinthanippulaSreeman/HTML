#include "udf.h"

/* Physical constants */
#define VOF_NUC    5e-4
#define R_B        1e-6
#define C_VAP      5.0
#define C_COND     0.01
#define P_OP       101325.0
#define T_INF      323.15
#define PV_REF     12352.0        // Vapor pressure at 50°C [Pa]
#define L_EV       2.382e6        // J/kg
#define CPL        4181.5         // J/kg·K
#define LAMBDA_L   0.64057        // W/m·K
#define R_V        461.5          // J/kg·K
#define DT         0.001          // s

DEFINE_LINEARIZED_MASS_TRANSFER(zgb_thermo_EZGB, cell, thread, from_index,
                                 from_species_index, to_index, to_species_index,
                                 d_mdot_d_vof_from, d_mdot_d_vof_to)
{
    Thread *liq = THREAD_SUB_THREAD(thread, from_index);
    Thread *vap = THREAD_SUB_THREAD(thread, to_index);

    real rho_l = C_R(cell, liq);       // Liquid density
    real rho_v = C_R(cell, vap);       // Vapor density
    real vof_v = C_VOF(cell, vap);     // VOF of vapor
    real T     = C_T(cell, thread);    // Temperature
    real press = C_P(cell, thread) + P_OP; // Pressure in cell + atmospheric pressure
    real k     = C_K(cell, thread);    // Thermal conductivity

    real m_dot = 0.0, dp = 0.0, m_source = 0.0;

    /* Basic sanity checks */
    if (rho_l <= 0.0 || rho_v <= 0.0 || fabs(rho_l - rho_v) < 1e-6)
        return 0.0;  // Return 0 if densities are not physically meaningful

    real r_rho_lv = 1.0 / rho_v - 1.0 / rho_l; // Reciprocal density difference
    real rho_m = rho_l * (1.0 - vof_v) + rho_v * vof_v;  // Mixture density
    real p_tur = 0.39 * k * rho_m;  // Turbulence term for pressure

    /* Vapor pressure calculation with temperature */
    real slope = (L_EV / T) * (rho_l * rho_v / (rho_l - rho_v));
    real Pv_T = PV_REF + slope * (T - T_INF) + 0.5 * p_tur;

    /* Thermal suppression factor */
    real alpha_l = LAMBDA_L / (rho_l * CPL);
    real sqrt_t = sqrt(DT);

    if (press < Pv_T)  // Evaporation condition
    {
        dp = MAX(Pv_T - press, 1e-4);  // Ensure dp isn't too small
        real suppression = (rho_l * CPL * sqrt(alpha_l) * (T_INF - T)) / (rho_v * L_EV * sqrt_t);
        real driving_term = sqrt((2.0 / 3.0) * dp / rho_l) - suppression;

        if (driving_term > 0.0)
        {
            m_dot = C_VAP * (3.0 * VOF_NUC * (1.0 - vof_v) * rho_v / R_B) * driving_term;
        }
    }
    else  // Condensation condition
    {
        dp = MAX(press - Pv_T, 1e-4);  // Ensure dp isn't too small
        real suppression = (rho_l * CPL * sqrt(alpha_l) * (T - T_INF)) / (rho_v * L_EV * sqrt_t);
        real driving_term = sqrt((2.0 / 3.0) * dp / rho_l) - suppression;

        if (driving_term > 0.0)
        {
            m_dot = -C_COND * (3.0 * vof_v * rho_v / R_B) * driving_term;
        }
    }

    /* Ensure mass flow rate is within physical bounds */
    if (isnan(m_dot) || isinf(m_dot) || fabs(m_dot) > 1e5)
        m_dot = 0.0;

    m_source = m_dot;
    *d_mdot_d_vof_from = m_dot;
    *d_mdot_d_vof_to   = -m_dot;

    /* Optional: Update storage if needed (e.g., for diagnostics) */
    if (NNULLP(THREAD_STORAGE(thread, SV_MT_DS_DP)))
    {
        C_STORAGE_R(cell, thread, SV_MT_DS_DP) = fabs(r_rho_lv * m_source / dp);
    }

    return m_source;  // Return the mass source rate
}
