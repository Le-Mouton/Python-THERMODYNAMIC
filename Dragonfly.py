import model

aero = model.Aerodynamique(6.2, 0.7, 0.03, -2.25, -10, 19.25, 1.5, 17.08, 6350 * 9.80665)
moteur = model.Moteur(1.005, 6.9, 1250.0, 1.4, 0.95, 0.8438, 0.8166, 0.9, 41376, 287.0528)
atmosphere = model.Atmosphere()

