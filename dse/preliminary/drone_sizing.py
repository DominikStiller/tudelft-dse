import numpy as np

rho_titanium = 4429
rho_cfrp = 1750
sigma_y_titanium = 900e6
sigma_ult_titanium = sigma_y_titanium / 1.5
sigma_ult_cfrp = 0.3 * 3e9


def m_to_ft(m):
    return m * 3.281


def ft2_to_m2(ft):
    return ft * 3.281**2


def kg_to_lbs(kg):
    return kg * 2.205


def lbs_to_kg(lbs):
    return lbs / 2.205


def calc_spar(rho, sigma_y, T_per_rotor, length):
    sigma_design = 0.9 * sigma_y
    t = 2e-3

    diameter_to_length_ratio = np.sqrt(4 * T_per_rotor / (np.pi * sigma_design * t * length))
    diameter = diameter_to_length_ratio * length
    volume = np.pi * t * diameter * length
    return volume * rho, diameter


def calc_blade_mass(rho, chord, radius):
    t = 0.05 * chord  # for clf5605, t/c=0.05
    A = t * chord / 2  # triangle
    V = A * radius
    return rho * V


def weights(
    n_blades,
    n_rotors,
    coaxial,
    n_legs,
    radius,
    c_to_R_ratio,
    tip_speed,
    MTOM,
    wet_area_fuselage,
    fuselage_length,
    engine_mass,
    spar_mass,
    blade_mass,
):
    """
    Calculate mass of a drone.

    Args:
        n_blades: number of blades per rotor
        n_rotors: number of rotors (if coaxial, number of axes)
        coaxial: whether the rotors are coaxial
        n_legs: number of legs
        radius: blade radius
        c_to_R_ratio: ratio of chord to radius
        tip_speed: tip speed
        MTOM: maximum take-off mass
        wet_area_fuselage: wetted area of the fuselage
        fuselage_length: length of the fuselage
        engine_mass: mass of a single engine
        spar_mass: mass of rotor spars
        blade_mass: mass of rotor blades

    Returns:
        Masses of components
    """
    n_engines = n_rotors
    n_total_rotors = 2 * n_rotors if coaxial else n_rotors

    # Convert to imperial
    c, R, OR = m_to_ft(np.array([radius * c_to_R_ratio, radius, tip_speed]))
    GW = kg_to_lbs(MTOM)
    Sw = ft2_to_m2(wet_area_fuselage)

    g_E = m_to_ft(9.81)
    J = c * R**3 / 3

    # Obtain weights in lbs
    W_rotors = n_total_rotors * n_blades * kg_to_lbs(blade_mass)
    W_engines = kg_to_lbs(n_engines * engine_mass)
    W_hubs = (
        n_rotors
        * (
            0.0037
            * n_blades**0.28
            * R**1.5
            * OR**0.43
            * (2 / 3 * W_rotors + g_E * J / R**2) ** 0.55
        )
        / 5
    )
    W_fuselage = 6.9 * (GW / 1000) ** 0.49 * fuselage_length**0.61 * Sw**0.25
    W_spars = kg_to_lbs(n_rotors * spar_mass)
    W_legs = 40 * (GW / 1000) ** (2 / 3) * n_legs**0.54
    W_misc = kg_to_lbs(200)  # 200 kg for miscellaneous items

    return lbs_to_kg(
        np.array(
            [
                W_rotors,
                W_engines,
                W_hubs,
                W_fuselage,
                W_spars,
                W_legs,
                W_misc,
            ]
        )
    )


def iterate_weights(
    payload_mass,
    fuel_mass,
    n_blades,
    n_rotors,
    coaxial,
    n_legs,
    radius,
    c_to_R_ratio,
    tip_speed,
    MTOM,
    wet_area_fuselage,
    fuselage_length,
    engine_mass,
    spar_mass,
    blade_mass,
):
    diff = 1
    while diff > 0.01:
        W_array = weights(
            n_blades,
            n_rotors,
            coaxial,
            n_legs,
            radius,
            c_to_R_ratio,
            tip_speed,
            MTOM,
            wet_area_fuselage,
            fuselage_length,
            engine_mass,
            spar_mass,
            blade_mass,
        )
        MTOM_new = np.sum(W_array) + payload_mass + fuel_mass
        diff = abs((MTOM_new - MTOM) / MTOM)
        MTOM = MTOM_new

    return MTOM, W_array


if __name__ == "__main__":
    # diff = 10
    # MTOM = 3000
    # useful_mass = 800  # Payload + fuel
    # R = 15
    # tip_speed = 200
    # Sw = 117
    # m_e = 600
    # while diff > 0.01:
    #     W_array = weights(
    #         n_blades=6,
    #         n_engines=2,
    #         n_legs=2,
    #         radius=R,
    #         tip_speed=tip_speed,
    #         MTOM=MTOM,
    #         wet_area_fuselage=Sw,
    #         engine_mass=m_e,
    #     )
    #     diff = abs(((np.sum(W_array) + useful_mass) - MTOM) / MTOM)
    #     MTOM = np.sum(W_array) + useful_mass
    #
    # print(W_array)
    # print(np.sum(W_array))

    # (n_rotors, n_blades, radius)
    params = [
        # (2, 6, 9.5),
        # (3, 6, 8),
        # (4, 3, 9),
        # (4, 4, 8.5),
        # (4, 6, 7),
        # (6, 6, 6),
        # (8, 3, 7),
        # (8, 4, 6),
        # (8, 6, 5.5),
        (4, 6, 7),
    ]
    for n_rotors, n_blades, radius in params:
        rho_spar = rho_cfrp
        sigma_ult_spar = sigma_ult_cfrp
        # rho_spar = rho_titanium
        # sigma_ult_spar = sigma_ult_titanium

        rho_blades = rho_cfrp

        payload_mass = 340
        fuel_mass = 500
        coaxial = True
        n_legs = 3
        c_to_R_ratio = 1 / 20
        tip_speed = 220 * 0.85
        MTOM = 2700
        wet_area_fuselage = 5 * 3 * 2
        fuselage_length = 5
        engine_mass = 40
        spar_length = radius * 1.2
        spar_mass, spar_diameter = calc_spar(
            rho_spar, sigma_ult_spar, 11130 / n_rotors, spar_length
        )
        blade_mass = calc_blade_mass(rho_blades, radius * c_to_R_ratio, radius)

        W_total, W_components = iterate_weights(
            payload_mass,
            fuel_mass,
            n_blades,
            n_rotors,
            coaxial,
            n_legs,
            radius,
            c_to_R_ratio,
            tip_speed,
            MTOM,
            wet_area_fuselage,
            fuselage_length,
            engine_mass,
            spar_mass,
            blade_mass,
        )
        print(f"{n_rotors=} {n_blades=} {spar_diameter=:.2} {W_total=}")
        print(W_components)
