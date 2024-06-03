import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar, minimize_scalar


class Moteur:

    def __init__(self, Cp, facteur_compression, temperature_max, gamma, rendement_tuyere, rendement_turbine,
                 rendement_combustion, rapport_compression, combustion_energie, r):

        self.Cp = Cp
        self.facteur_compression = facteur_compression
        self.temperature_max = temperature_max
        self.gamma = gamma

        self.rendement_turbine = rendement_turbine
        self.rendement_tuyere = rendement_tuyere
        self.rendement_combustion = rendement_combustion
        self.rapport_compression = rapport_compression
        self.combustion_energie = combustion_energie
        self.r = r

    def setCondition(self, Pext, Text, Vavion):

        self.Vavion = Vavion
        self.Point1 = {"enthalpie": self.Cp * Text, "température": Text, "pression": Pext}
        return self.Point1

    def point2(self):

        H2 = self.Point1["enthalpie"] + (self.Vavion ** 2) / 2 / 1000
        T2 = H2 / self.Cp
        P2 = self.Point1["pression"] * (T2 / self.Point1["température"]) ** (self.gamma / (self.gamma - 1))
        self.Point2 = {"enthalpie": H2, "température": T2, "pression": P2}
        return self.Point2

    def point3(self):

        P3 = self.Point2["pression"] * self.facteur_compression
        T3 = self.Point2["température"] * (P3 / self.Point2["pression"]) ** ((self.gamma - 1) / self.gamma)
        T3_reel = self.Point2["température"] + (T3 - self.Point2["température"]) / self.rendement_combustion
        H3 = T3_reel * self.Cp
        self.Point3 = {"enthalpie": H3, "température": T3_reel, "pression": P3}
        return self.Point3

    def point4(self):

        T4 = self.temperature_max
        P4 = self.Point3["pression"] * self.rapport_compression
        H4 = T4 * self.Cp
        self.Point4 = {"enthalpie": H4, "température": T4, "pression": P4}
        return self.Point4

    def point5(self):

        H5 = self.Point4["enthalpie"] - (self.Point3["enthalpie"] - self.Point2["enthalpie"])
        T5 = H5 / self.Cp
        T5_reel = self.Point4["température"] + (T5 - self.Point4["température"]) / self.rendement_turbine
        P5 = self.Point4["pression"] * (T5_reel / self.Point4["température"]) ** (self.gamma / (self.gamma - 1))
        self.Point5 = {"enthalpie": H5, "température": T5, "pression": P5}
        return self.Point5

    def point6(self):

        T6 = self.Point5["température"] * (self.Point1["pression"] / self.Point5["pression"]) ** (
                (self.gamma - 1) / self.gamma)
        T6_reel = self.Point5["température"] + (T6 - self.Point5["température"]) * self.rendement_tuyere
        H6 = T6_reel * self.Cp
        P6 = self.Point1["pression"]
        self.Point6 = {"enthalpie": H6, "température": T6_reel, "pression": P6}
        return self.Point6

    def points(self):

        print(self.Point1, "\n",
              self.point2(), "\n",
              self.point3(), "\n",
              self.point4(), "\n",
              self.point5(), "\n",
              self.point6(), "\n")

    def performance(self):

        # Points calculés
        self.point2()
        self.point3()
        self.point4()
        self.point5()
        self.point6()

        qc = self.Point4["enthalpie"] - self.Point3["enthalpie"]
        fuel_air_ratio = qc / self.combustion_energie

        vitesse_ejection = (2 * (self.Point5["enthalpie"] - self.Point6["enthalpie"]) * 1000) ** 0.5

        col_diametre = 29.13e-2
        section_col = (col_diametre ** 2 * 3.14159) / 4
        debit_air = 0.0404 * section_col * self.Point5["pression"] / ((self.Point5['température']) ** 0.5) * 1000

        thrust = debit_air * (vitesse_ejection - self.Vavion) / 1000

        debit_fuel = debit_air * fuel_air_ratio * 1000

        conso_spe = debit_fuel / thrust

        rendement_moteur = thrust * self.Vavion / (debit_fuel * self.combustion_energie * 1000)

        return {
            "vitesse ejection": vitesse_ejection,
            "debit air": debit_air,
            "debit fuel": debit_fuel,
            "ratio": fuel_air_ratio,
            "poussee": thrust,
            "consomation specifique": conso_spe,
            "rendement moteur": rendement_moteur
        }

    def plot_diagramme_PT(self):

        points = [
            self.Point1,
            self.point2(),
            self.point3(),
            self.point4(),
            self.point5(),
            self.point6(),
            self.Point1
        ]

        temperatures = [point["température"] for point in points]
        pressions = [point["pression"] for point in points]

        plt.figure(figsize=(10, 6))
        plt.plot(temperatures, pressions, 'o-', label='Cycle thermodynamique', color='blue')
        plt.xlabel('Température (K)')
        plt.ylabel('Pression (kPa)')
        plt.title('Diagramme P-T (Pression-Température)')
        plt.grid(True)

        for i in range(len(points) - 1):
            T1, P1 = points[i]["température"], points[i]["pression"]
            T2, P2 = points[i + 1]["température"], points[i + 1]["pression"]
            plt.plot([T1, T2], [P1, P2], 'r')

        n_lines = 10
        P = np.linspace(min(pressions), max(pressions), 50)

        T_initial = np.linspace(min(temperatures), max(temperatures), n_lines)

        for T0 in T_initial:
            adiabat = T0 * (P / P[0]) ** ((self.gamma - 1) / self.gamma)
            plt.plot(adiabat, P, 'g--')

        plt.xscale("linear")
        plt.yscale("linear")
        plt.xlim(min(temperatures) - 50, max(temperatures) + 50)
        plt.ylim(min(pressions) - 10, max(pressions) + 10)

        plt.legend()
        plt.show()


class Atmosphere:
    def __init__(self):
        self.T0 = 288.15
        self.P0 = 101.325
        self.g = 9.80665
        self.R = 287.0528
        self.alpha = 6.5 / 1000
        self.TropopauseAltitude = 11000
        self.TropopauseTemperature = 216.65
        self.TropopausePression = 22.6377
        self.gamma = 1.4

    def get_conditions(self, altitude):

        if altitude <= self.TropopauseAltitude:
            temperature = self.T0 - self.alpha * altitude
            pression = self.P0 * (1 - (self.alpha * altitude) / self.T0) ** (self.g / (self.alpha * self.R))
        else:
            temperature = self.TropopauseTemperature
            pression = self.TropopausePression * np.exp(
                -self.g * (altitude - self.TropopauseAltitude) / (self.R * temperature))

        masse_vol = pression * 1000 / (self.R * temperature)
        vitesse_son = (self.gamma * self.R * temperature) ** 0.5

        return {
            "pression": pression,
            "temperature": temperature,
            "masse_vol": masse_vol,
            "vitesse_son": vitesse_son
        }

    def plot_p_x_altitude(self, altitudes):

        pressions = [self.get_conditions(alt)["pression"] for alt in altitudes]

        plt.figure(figsize=(8, 6))
        plt.plot(pressions, altitudes, 'b-', label='Pression x Altitude')
        plt.ylabel('Altitude (m)')
        plt.xlabel('Pression (kPa)')
        plt.title('Diagramme Pression-Altitude')
        plt.grid(True)
        plt.legend()
        plt.show()

    def plot_t_x_altitude(self, altitudes):

        temperatures = [self.get_conditions(alt)["temperature"] for alt in altitudes]

        plt.figure(figsize=(8, 6))
        plt.plot(temperatures, altitudes, 'r-', label='Température x Altitude')
        plt.ylabel('Altitude (m)')
        plt.xlabel('Température (K)')
        plt.title('Diagramme Température-Altitude')
        plt.grid(True)
        plt.legend()
        plt.show()


class Aerodynamique:
    def __init__(self, AR, e_oswald, Cd0, i0, imin, idec, Cl_max, surface_alaire, poids_avion):

        self.AR = AR
        self.e_oswald = e_oswald
        self.Cd0 = Cd0
        self.i0 = i0
        self.imin = imin
        self.idec = idec
        self.Cl_max = Cl_max
        self.surface_alaire = surface_alaire
        self.poids_avion = poids_avion

    def coefficient_portance(self, incidence):

        if self.imin <= incidence <= self.idec:
            Cl = self.Cl_max * np.sin((np.pi / 2) * (incidence - self.imin) / (self.idec - self.imin)) * self.AR / (
                        self.AR + 2)
            return Cl
        else:
            return self.Cl_max

    def coefficient_trainee(self, incidence):

        Cl = self.coefficient_portance(incidence)
        Cd = self.Cd0 + Cl ** 2 / (np.pi * self.AR * self.e_oswald)
        return Cd

    def finesse(self, incidence):

        Cl = self.coefficient_portance(incidence)
        Cd = self.coefficient_trainee(incidence)
        finesse = Cl / Cd
        return finesse

    def incidence_pour_portance(self, Cl):

        incidence = Cl / (2 * np.pi) * 180 / np.pi + self.i0
        return incidence

    def incidence_pour_trainee(self, Cd):

        incidences = [i for i in range(self.imin, int(self.idec) + 1)]
        for i in incidences:
            if np.isclose(self.coefficient_trainee(i), Cd, atol=1e-5):
                return i
        return None

    def incidence_pour_finesse(self, finesse):

        incidences = [i for i in range(self.imin, int(self.idec) + 1)]
        for i in incidences:
            if np.isclose(self.finesse(i), finesse, atol=1e-5):
                return i
        return None

    def finesse_max(self):

        incidences = np.linspace(self.imin, self.idec, 500)

        cls = [self.coefficient_portance(inc) for inc in incidences]
        cds = [self.coefficient_trainee(inc) for inc in incidences]
        finesse = [self.finesse(inc) for inc in incidences]

        max_finesse_index = np.argmax(finesse)
        max_finesse = finesse[max_finesse_index]
        max_finesse_incidence = incidences[max_finesse_index]

        return {
            "incidence": max_finesse_incidence,
            "finesse_max": max_finesse,
            "cd_finesse_max": cds[max_finesse_index],
            "cl_finesse_max": cls[max_finesse_index]
        }

    def plot_cd_cl_x_incidence(self):

        incidences = np.linspace(self.imin, self.idec, 300)

        cls = [self.coefficient_portance(inc) for inc in incidences]
        cds = [self.coefficient_trainee(inc) for inc in incidences]
        clxcd = [np.sqrt(cl) / cd for cl, cd in zip(cls, cds)]
        finesse = [self.finesse(inc) for inc in incidences]

        param = self.finesse_max()

        plt.figure(figsize=(12, 9))
        plt.plot(incidences, cls, label='Coefficient de portance (Cl)', color='blue')
        plt.plot(incidences, clxcd, label='Rapport √(Cl)/Cd', color='purple')
        plt.plot(incidences, cds, label='Coefficient de traînée (Cd)', color='red')
        plt.plot(incidences, finesse, label='Finesse', color='green')

        plt.plot(param["incidence"], param["finesse_max"], 'o', label='Finesse max', color='black', markersize=10)

        plt.title('Coefficients de portance, de traînée, rapport √(Cl)/Cd et finesse en fonction de l\'Incidence')
        plt.xlabel('Angle d\'incidence (degrés)')
        plt.ylabel('Coefficient')
        plt.legend()
        plt.grid(True)
        plt.savefig("all_incidence.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(incidences, cls, color='blue')
        plt.title('Coefficients de Portance en fonction de l\'Incidence')
        plt.xlabel('Angle d\'incidence (degrés)')
        plt.ylabel('Coefficient de portance (Cl)')
        plt.legend()
        plt.grid(True)
        plt.savefig("cl_incidence.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(incidences, cds, color='red')
        plt.title('Coefficients de Trainée en fonction de l\'Incidence')
        plt.xlabel('Angle d\'incidence (degrés)')
        plt.ylabel('Coefficient de trainée (Cl)')
        plt.legend()
        plt.grid(True)
        plt.savefig("cd_incidence.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(incidences, clxcd, color='purple')
        plt.title('Rapport √(Cl)/Cd en fonction de l\'Incidence')
        plt.xlabel('Angle d\'incidence (degrés)')
        plt.ylabel('Rapport √(Cl)/Cd')
        plt.legend()
        plt.grid(True)
        plt.savefig("clxcd_incidence.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(incidences, finesse, color='green')
        plt.title('Finesse en fonction de l\'Incidence')
        plt.xlabel('Angle d\'incidence (degrés)')
        plt.ylabel('Finesse')
        plt.legend()
        plt.grid(True)
        plt.savefig("finesse_incidence.png")
        plt.show()


class Performance:
    def __init__(self, moteur, atmosphere, aerodynamique):
        self.moteur = moteur
        self.atmosphere = atmosphere
        self.aerodynamique = aerodynamique

    def calc_vitesse(self, altitude):
        conditions_atm = self.atmosphere.get_conditions(altitude)
        fmax = self.aerodynamique.finesse_max()

        def f_vitesse(V):
            Mach = V / conditions_atm['vitesse_son']
            portance = (0.5 * conditions_atm['masse_vol'] * V ** 2 * self.aerodynamique.surface_alaire * fmax[
                "cl_finesse_max"]) / np.sqrt(1 - Mach ** 2)
            return portance - self.aerodynamique.poids_avion

        result = root_scalar(f_vitesse, bracket=(0, 0.9 * conditions_atm['vitesse_son']), method='bisect')
        vitesse = result.root
        return vitesse

    def plot_vitesse_mach_vs_altitude(self, max_altitude=14000, step=50):
        altitudes = np.arange(2000, max_altitude, step)
        vitesses = []
        mach_numbers = []
        for altitude in altitudes:
            conditions = self.atmosphere.get_conditions(altitude)
            vitesse = self.calc_vitesse(altitude)
            mach_number = vitesse / conditions['vitesse_son']
            vitesses.append(vitesse)
            mach_numbers.append(mach_number)

        plt.figure(figsize=(12, 9))

        plt.plot(vitesses, altitudes / 1000, 'b-')
        plt.title("Vitesse en fonction de l'altitude")
        plt.xlabel('Vitesse (m/s)')
        plt.ylabel('Altitude (km)')
        plt.grid(True)
        plt.savefig("Vitesse en fonction de l'altitude.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(mach_numbers, altitudes / 1000, 'r-')
        plt.title("Nombre de Mach en fonction de l'altitude")
        plt.xlabel('Nombre de Mach')
        plt.ylabel('Altitude (km)')
        plt.grid(True)

        plt.tight_layout()
        plt.savefig("Mach en fonction de l'altitude.png")
        plt.show()

    def plot_poussemax_rendement(self):
        condition_atm = self.atmosphere.get_conditions(12000)
        vitesses = np.linspace(0, 0.7 * condition_atm['vitesse_son'], 500)
        poussée_max = []
        rendement_propulsif = []

        for v in vitesses:
            self.moteur.setCondition(condition_atm['pression'], condition_atm['temperature'], v)
            perf_moteur = self.moteur.performance()
            poussée_max.append(perf_moteur['poussee'])
            rendement_propulsif.append(perf_moteur['rendement moteur'] * 10e7)

        plt.figure(figsize=(12, 9))
        plt.plot(vitesses, poussée_max, 'b-')
        plt.title("Évolution de la poussée maximale en fonction de la vitesse de l'Avion")
        plt.xlabel('Vitesse de l\'avion (m/s)')
        plt.ylabel('Poussée maximale (kN)')
        plt.grid(True)
        plt.savefig("poussemax_vitesse.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(vitesses, rendement_propulsif, 'r-')
        plt.title("Évolution du rendement moteur en fonction de la vitesse de l'Avion")
        plt.xlabel('Vitesse de l\'avion (m/s)')
        plt.ylabel('Rendement (%)')
        plt.grid(True)
        plt.savefig("rendement_vitesse.png")
        plt.show()

    def plot_performance_vs_altitude(self):
        altitudes = np.linspace(2000, 14000, 500)
        poussée_max = []
        debit_fuel = []
        conso_specifique = []

        for altitude in altitudes:
            condition_atm = self.atmosphere.get_conditions(altitude)
            vitesse = 0.7 * condition_atm['vitesse_son']
            self.moteur.setCondition(condition_atm['pression'], condition_atm['temperature'], vitesse)
            perf_moteur = self.moteur.performance()

            poussée_max.append(perf_moteur['poussee'])
            debit_fuel.append(perf_moteur['debit fuel'])
            conso_specifique.append(perf_moteur['consomation specifique'])

        plt.figure(figsize=(12, 9))
        plt.plot(poussée_max, altitudes / 1000, 'b-')
        plt.title("Évolution de la poussée maximale en fonction de l'altitude")
        plt.ylabel('Altitude (km)')
        plt.xlabel('Poussée maximale (kN)')
        plt.grid(True)
        plt.savefig("pousseemaximale_altitude.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(poussée_max, debit_fuel, 'g-')
        plt.title("Évolution du débit de fuel en fonction de la poussée maximale")
        plt.xlabel('Poussée maximale (kN)')
        plt.ylabel('Débit de fuel (g/s)')
        plt.grid(True)
        plt.savefig("debitfuel_poussee.png")
        plt.show()

        plt.figure(figsize=(12, 9))
        plt.plot(poussée_max, conso_specifique, 'r-')
        plt.title("Évolution de la consommation spécifique en fonction de la poussée maximale")
        plt.xlabel('Poussée maximale (kN)')
        plt.ylabel('Consommation spécifique (g/s/kN)')
        plt.grid(True)
        plt.savefig("consospe_pousseemaximale.png")
        plt.show()

    def calc_trainee_minimale(self, altitude):
        conditions_atm = self.atmosphere.get_conditions(altitude)
        masse_vol = conditions_atm['masse_vol']
        fmax = self.aerodynamique.finesse_max()
        cl_finesse_max = fmax["cl_finesse_max"]
        cd_finesse_max = fmax["cd_finesse_max"]

        vitesse = (2 * self.aerodynamique.poids_avion / (
                masse_vol * self.aerodynamique.surface_alaire * cl_finesse_max)) ** 0.5
        Mach = vitesse / conditions_atm['vitesse_son']
        cl_corrige = cl_finesse_max / (1 + 0.2 * Mach ** 2)  # Correction de Cl par le Mach
        cd_corrige = cd_finesse_max * (1 + 0.2 * Mach ** 2)  # Correction de Cd par le Mach

        trainee = 0.5 * masse_vol * vitesse ** 2 * self.aerodynamique.surface_alaire * cd_corrige
        return trainee

    def calc_poussee_maximale(self, altitude):
        conditions_atm = self.atmosphere.get_conditions(altitude)
        vitesse = 0
        self.moteur.setCondition(conditions_atm['pression'], conditions_atm['temperature'], vitesse)
        perf_moteur = self.moteur.performance()
        return 2 * perf_moteur['poussee'] * 1000

    def plafond(self):
        def difference_poussee_trainee(altitude):
            trainee_minimale = self.calc_trainee_minimale(altitude)
            poussee_maximale = self.calc_poussee_maximale(altitude)
            return abs(poussee_maximale - trainee_minimale)

        res = minimize_scalar(difference_poussee_trainee, bounds=(0, 20000), method='bounded')
        return res.x, difference_poussee_trainee(res.x)

    def calc_trainee(self, altitude, vitesse):
        conditions_atm = self.atmosphere.get_conditions(altitude)
        masse_vol = conditions_atm['masse_vol']

        cl = (2 * self.aerodynamique.poids_avion) / (masse_vol * self.aerodynamique.surface_alaire * vitesse ** 2)

        cd = self.aerodynamique.Cd0 + cl ** 2 * 1 / (np.pi * self.aerodynamique.e_oswald * self.aerodynamique.AR)

        Mach = vitesse / conditions_atm['vitesse_son']
        cd_corrige = cd * (1 + 0.2 * Mach ** 2)

        trainee = 0.5 * masse_vol * vitesse ** 2 * self.aerodynamique.surface_alaire * cd_corrige
        return trainee / 1000

    def calc_poussee_max(self, vitesse, altitude):
        conditions_atm = self.atmosphere.get_conditions(altitude)
        self.moteur.setCondition(conditions_atm['pression'], conditions_atm['temperature'], vitesse)
        perf_moteur = self.moteur.performance()
        return 2 * perf_moteur['poussee']

    def plot_trainee_pousse(self, altitude):
        vitesses = np.linspace(130, 0.9 * self.atmosphere.get_conditions(altitude)['vitesse_son'], 500)
        trainees = []
        poussees = []

        for vitesse in vitesses:
            trainee_list = []
            pousse = self.calc_poussee_max(vitesse, altitude)
            trainee_list.append(self.calc_trainee(altitude, vitesse))
            trainees.append(max(trainee_list))  # Prendre la valeur maximale de la traînée pour chaque vitesse
            poussees.append(pousse)

        plt.figure(figsize=(12, 9))
        plt.plot(vitesses * 3.6, trainees, 'k-', label='Traînée (kN)')  # Conversion de la vitesse en km/h
        plt.plot(vitesses * 3.6, poussees, 'r-',
                 label='Poussée max des moteurs (kN)')  # Conversion de la vitesse en km/h
        plt.title(f'z = {altitude / 1000} km')
        plt.xlabel('V (km/h)')
        plt.ylabel('T (kN)')
        plt.legend()
        plt.xlim(450, 950)
        plt.ylim(0, 10)
        plt.grid(True)
        plt.show()
