import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Show diagram?
DiagrammKraftverlauf = 1
DiagrammLeerlaufkraft = 1

DateiAuswertegrenzen = 'Auswertegrenzen_Kraft_Zerspanungsversuche.txt'
Auswertegrenzen = np.loadtxt(DateiAuswertegrenzen)


# Initialisation averaged process forces
Versuche = len(Auswertegrenzen) * [0]
Schnittkraft = len(Auswertegrenzen) * [0.0]
Passivkraft = len(Auswertegrenzen) * [0.0]
Vorschubkraft = len(Auswertegrenzen) * [0.0]
# Standardabweichungen
SigmaSchnittkraft = len(Auswertegrenzen) * [0.0]
SigmaPassivkraft = len(Auswertegrenzen) * [0.0]
SigmaVorschubkraft = len(Auswertegrenzen) * [0.0]

Pfad = r"D:\70_Publications\2021\03_Hobeln\50_Experimente\01_Kraftmessungen\Gruppe1"  # Windows
# Pfad =
# r"/Users/fabiankneubuhler/Desktop/Kraftrauswertung/Python/Daten/Gruppe1"
# # Mac

Dateien = []


# r=root, d=directories, f = files
# (https://www.mkyong.com/python/python-how-to-list-all-files-in-a-directory/)
for r, d, f in os.walk(Pfad):
    for Datei in f:
        if ('.dat' in Datei) and ('ProzKraft' in Datei):
            Dateien.append(os.path.join(r, Datei))


# Kräfte ohne Eingriff (Eigenmasse der MeÃŸplattform usw.)
Fc_Leerlauf_vorne = len(Auswertegrenzen) * [0]  # aus vorderem Teil ermitteln
Fc_Leerlauf_hinten = len(Auswertegrenzen) * [0]  # aus hinterem Teil ermitteln
Fp_Leerlauf_vorne = len(Auswertegrenzen) * [0]
Fp_Leerlauf_hinten = len(Auswertegrenzen) * [0]
Ff_Leerlauf_vorne = len(Auswertegrenzen) * [0]
Ff_Leerlauf_hinten = len(Auswertegrenzen) * [0]
#
FfLeerAvg = len(Auswertegrenzen) * [np.nan]
FcLeerAvg = len(Auswertegrenzen) * [np.nan]
FfLeerStd = len(Auswertegrenzen) * [np.nan]
FcLeerStd = len(Auswertegrenzen) * [np.nan]


for Datei in Dateien:

    # Read head of "ProzKraft.dat" file (number of experiment, ...)
    D2 = open(Datei, "r")
    Zeile = D2.readline()
    D2.close()

    # Find number of experiment
    VStart = Zeile.find(':') + 1
    VEnde = Zeile.find(',')
    Versuchsnummer = Zeile[VStart:VEnde]

    # Find cutting speed
    #vcStart = Zeile.find(':', VStart)+1
    #vcEnde = Zeile.find(',', VEnde+1)
    #vc = Zeile[vcStart:vcEnde]

    # Find feed
    #fStart = Zeile.find(':', vcStart)+1
    #fEnde = Zeile.find(',', vcEnde+1)
    #f = Zeile[fStart:fEnde]
    #
    Daten = np.loadtxt(Datei, skiprows=2)
    Fc = Daten[:, 1]  # Cutting force
    Ff = Daten[:, 3]  # Feed force
    Fp = Daten[:, 2]  # Passive force

    AnzWerte = Daten[:, 1].size

    #BeginnMittelung = int(0.2*AnzWerte)
    #EndeMittelung = int(0.8*AnzWerte)

    BeginnMittelung = int(
        Auswertegrenzen[int(Versuchsnummer) - 1, 1] / 100 * AnzWerte)
    EndeMittelung = int(
        Auswertegrenzen[int(Versuchsnummer) - 1, 2] / 100 * (AnzWerte - 1))

    Wandstaerke = float(Auswertegrenzen[int(Versuchsnummer) - 1, 5])
    #
    # AvgSchnittkraft = np.mean(Daten[BeginnMittelung:EndeMittelung,2])/Wandstaerke # Falsch! Koordinatensytem der KraftmeÃŸplattform ist nicht gleich! HK, Fr, 05.06.2020
    # AvgPassivkraft = np.mean(Daten[BeginnMittelung:EndeMittelung,1])/Wandstaerke # x-Achse ist auf KraftmeÃŸplattform und Maschine gleich! HK, Fr, 05.06.2020
    # AvgVorschubkraft = np.mean(Daten[BeginnMittelung:EndeMittelung,3])/Wandstaerke # Falsch! Koordinatensytem der KraftmeÃŸplattform ist nicht gleich! HK, Fr, 05.06.2020
    # AvgVorschubkraft = np.mean(Daten[BeginnMittelung:EndeMittelung,2])/Wandstaerke # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
    # AvgSchnittkraft =
    # np.mean(Daten[BeginnMittelung:EndeMittelung,3])/Wandstaerke #
    # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
    # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
    AvgVorschubkraft = np.mean(Ff[BeginnMittelung:EndeMittelung]) / Wandstaerke
    # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
    AvgSchnittkraft = np.mean(Fc[BeginnMittelung:EndeMittelung]) / Wandstaerke
    #########
    #
    Fc_Leerlauf_vorne[int(
        Versuchsnummer) - 1] = np.mean(Daten[round(0.999 * len(Daten[:, 0])):-1, 3])
    Fc_Leerlauf_hinten[int(Versuchsnummer) - 1] = np.mean(
        Daten[0:round(0.001 * len(Daten[:, 0])), 3])  # aus hinterem Teil ermitteln
    #Fp_Leerlauf_vorne = []
    #Fp_Leerlauf_hinten = []
    Ff_Leerlauf_vorne[int(
        Versuchsnummer) - 1] = np.mean(Daten[round(0.999 * len(Daten[:, 0])):-1, 2])
    Ff_Leerlauf_hinten[int(
        Versuchsnummer) - 1] = np.mean(Daten[0:round(0.001 * len(Daten[:, 0])), 2])
    #
    # Varianz berechnen
    # SigSchnittkraft = np.sqrt(np.var(Daten[BeginnMittelung:EndeMittelung,2]))/Wandstaerke # Falsch! Koordinatensytem der KraftmeÃŸplattform ist nicht gleich! HK, Fr, 05.06.2020
    # SigPassivkraft = np.sqrt(np.var(Daten[BeginnMittelung:EndeMittelung,3]))/Wandstaerke # x-Achse ist auf KraftmeÃŸplattform und Maschine gleich! HK, Fr, 05.06.2020
    # SigVorschubkraft = np.sqrt(np.var(Daten[BeginnMittelung:EndeMittelung,3]))/Wandstaerke # Falsch! Koordinatensytem der KraftmeÃŸplattform ist nicht gleich! HK, Fr, 05.06.2020
    # Standardabweichung berechnen
    # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
    SigSchnittkraft = (np.std(Fc[BeginnMittelung:EndeMittelung])) / Wandstaerke
    # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
    SigVorschubkraft = (
        np.std(Ff[BeginnMittelung:EndeMittelung])) / Wandstaerke
    #
    # Aktuelle Schnittkraft speichern; spÃ¤ter in Datei ausschreiben
    Schnittkraft[int(Versuchsnummer) - 1] = AvgSchnittkraft
    #Passivkraft[int(Versuchsnummer)-1] = AvgPassivkraft
    Vorschubkraft[int(Versuchsnummer) - 1] = AvgVorschubkraft
    Versuche[int(Versuchsnummer) - 1] = int(Versuchsnummer)
    SigmaSchnittkraft[int(Versuchsnummer) - 1] = SigSchnittkraft
    #SigmaPassivkraft[int(Versuchsnummer)-1] = SigPassivkraft
    SigmaVorschubkraft[int(Versuchsnummer) - 1] = SigVorschubkraft
    #
    #MedianSchnittkraft = np.median(Daten[:,2])
    #MedianPassivkraft = np.median(Daten[:,3])
    # https://matplotlib.org/users/usetex.html
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

#    plt.plot(Daten[:,0], Daten[:,3]/Wandstaerke, color="blue", linewidth=1.0, linestyle="-", label=r'Schnittkraft [N/mm]') # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
#    #plt.plot(Daten[:,0], Daten[:,1]/Wandstaerke, color="red", linewidth=1.0, linestyle="-", label=r'Passivkraft [N/mm]') # x-Achse ist auf KraftmeÃŸplattform und Maschine gleich! HK, Fr, 05.06.2020
#    plt.plot(Daten[:,0], Daten[:,2]/Wandstaerke, color="red", linewidth=1.0, linestyle="-", label=r'Vorschubkraft [N/mm]') # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
#
    LinieAvgFc = [AvgSchnittkraft, AvgSchnittkraft]
#    #LinieAvgFp = [AvgPassivkraft, AvgPassivkraft]
    LinieAvgFf = [AvgVorschubkraft, AvgVorschubkraft]
#    plt.xlim(min(Daten[:,0]), max(Daten[:,0]))  # Grenzen (min/max) X-Achse
#    plt.plot([Daten[BeginnMittelung,0], Daten[EndeMittelung,0]], LinieAvgFc, color="darkblue", linewidth=0.5, linestyle="--", label=r'Schnittkraft\_Mittel: '+str(round(AvgSchnittkraft,1)) + '[N/mm], Sigma: ' + str(round(SigSchnittkraft,1)))
#    #plt.plot([Daten[BeginnMittelung,0], Daten[EndeMittelung,0]], LinieAvgFp, color="orange", linewidth=0.5, linestyle="--", label=r'Passivkraft\_Mittel: '+str(round(AvgPassivkraft,1)) + '[N/mm], Sigma: ' + str(round(SigPassivkraft,1)))
#    plt.plot([Daten[BeginnMittelung,0], Daten[EndeMittelung,0]], LinieAvgFf, color="orange", linewidth=0.5, linestyle="--", label=r'Vorschubkraft\_Mittel: '+str(round(AvgVorschubkraft,1)) + '[N/mm], Sigma: ' + str(round(SigVorschubkraft,1)))
#    plt.legend(loc="lower center")
#
#    plt.xlabel(r'Zeit [$s$]')
#    plt.ylabel(r'Kraft [N/mm]')
#    plt.title(r'Proze{\ss}kraftverlauf, Versuch ' + str(Versuchsnummer) + ', vc= ' + vc + ', f= ' + f)

    if (DiagrammKraftverlauf != 0):
        fig, ax = plt.subplots(constrained_layout=True)
        ax.plot(Daten[:, 0], Fc[:] / Wandstaerke, color="blue", linewidth=1.0, linestyle="-",
                label=r'Cutting force [N/mm]')  # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
        ax.plot(Daten[:, 0], Ff[:] / Wandstaerke, color="red", linewidth=1.0, linestyle="-",
                label=r'Feed force [N/mm]')  # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
        ax.plot(Daten[:, 0], Fp[:] / Wandstaerke, color="green",
                linewidth=1.0, linestyle="-", label=r'Passive force [N/mm]')

        ax.plot([Daten[BeginnMittelung, 0], Daten[EndeMittelung, 0]], LinieAvgFc, color="darkblue", linewidth=0.5, linestyle="--",
                label=r'Average cutting force: ' + str(round(AvgSchnittkraft, 1)) + '[N/mm], Sigma: ' + str(round(SigSchnittkraft, 1)))
        #ax.plot([Daten[BeginnMittelung,0], Daten[EndeMittelung,0]], LinieAvgFp, color="orange", linewidth=0.5, linestyle="--", label=r'Passivkraft\_Mittel: '+str(round(AvgPassivkraft,1)) + '[N/mm], Sigma: ' + str(round(SigPassivkraft,1)))
        ax.plot([Daten[BeginnMittelung, 0], Daten[EndeMittelung, 0]], LinieAvgFf, color="orange", linewidth=0.5, linestyle="--",
                label=r'Average feed force: ' + str(round(AvgVorschubkraft, 1)) + '[N/mm], Sigma: ' + str(round(SigVorschubkraft, 1)))

        ax.set_xlabel(r'Time [$s$]')
        ax.set_ylabel(r'Force [N/mm]')
        #ax.set_title(r'Proze{\ss}kraftverlauf, Versuch ' + str(Versuchsnummer) + ', vc= ' + vc + ', f= ' + f)
        ax.set_title(r'Process force, Experiment ' + str(Versuchsnummer))
        ax.legend(loc="lower center")
        # ax.xlim(min(Daten[:,0]), max(Daten[:,0]))  # Grenzen (min/max)
        # X-Achse
        ax.xaxis.limit_range_for_scale(0, 1)
        # Actual plotting code omitted
        ax.xaxis.set_major_formatter(mtick.PercentFormatter(max(Daten[:, 0])))

        def forward(x):
            return x

        def inverse(x):
            return x
        secax = ax.secondary_xaxis('top', functions=(forward, inverse))
        secax.set_xlabel(r'Time [$s$]')
        #
        plt.savefig("Prozeszkraefte_" + Versuchsnummer + ".png", dpi=300)
        # plt.savefig('Prozeszkraefte')
        #
        plt.show()
        #
        plt.close()
        #

        # Show plot of force evaluation area
        ZoomBeginn = Auswertegrenzen[int(
            Versuchsnummer) - 1, 3] / 100.0 * max(Daten[:, 0])
        ZoomEnde = Auswertegrenzen[int(
            Versuchsnummer) - 1, 4] / 100.0 * max(Daten[:, 0])

        plt.xlim(ZoomBeginn, ZoomEnde)  # Grenzen (min/max) X-Achse
        plt.plot(Daten[:, 0], Fc[:] / Wandstaerke, color="blue", linewidth=1.0, linestyle="-",
                 label=r'Cutting force [N/mm]')  # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
        # plt.plot(Daten[:,0], Daten[:,1]/Wandstaerke, color="red",
        # linewidth=1.0, linestyle="-", label=r'Passivkraft [N/mm]') # x-Achse
        # ist auf KraftmeÃŸplattform und Maschine gleich! HK, Fr, 05.06.2020
        plt.plot(Daten[:, 0], Ff[:] / Wandstaerke, color="red", linewidth=1.0, linestyle="-",
                 label=r'Feed force [N/mm]')  # Korrigierte Zuordnung der ProzeÃŸkrÃ¤fte, HK, Fr, 05.06.2020
        plt.plot(Daten[:, 0], Fp[:] / Wandstaerke, color="green",
                 linewidth=1.0, linestyle="-", label=r'Passive force [N/mm]')

        LinieAvgFc = [AvgSchnittkraft, AvgSchnittkraft]
        #LinieAvgFp = [AvgPassivkraft, AvgPassivkraft]
        LinieAvgFp = [AvgVorschubkraft, AvgVorschubkraft]
        plt.plot([Daten[BeginnMittelung, 0], Daten[EndeMittelung, 0]], LinieAvgFc, color="darkblue", linewidth=0.5, linestyle="--",
                 label=r'Average cutting force: ' + str(round(AvgSchnittkraft, 1)) + '[N/mm], Sigma: ' + str(round(SigSchnittkraft, 1)))
        #plt.plot([Daten[BeginnMittelung,0], Daten[EndeMittelung,0]], LinieAvgFp, color="orange", linewidth=0.5, linestyle="--", label=r'Passivkraft\_Mittel: '+str(round(AvgPassivkraft)) + '[N/mm], Sigma: ' + str(round(SigPassivkraft)))
        plt.plot([Daten[BeginnMittelung, 0], Daten[EndeMittelung, 0]], LinieAvgFp, color="orange", linewidth=0.5, linestyle="--",
                 label=r'Average feed force: ' + str(round(AvgVorschubkraft, 1)) + '[N/mm], Sigma: ' + str(round(SigVorschubkraft, 1)))
        plt.legend(loc="lower center")

        plt.xlabel(r'Time [$s$]')
        plt.ylabel(r'Force [N/mm]')
        #plt.title(r'Proze{\ss}kraftverlauf, Versuch ' + str(Versuchsnummer) + ', vc= ' + vc + ', f= ' + f)
        plt.title(r'Process force, Experiment ' + str(Versuchsnummer))
        #
        plt.savefig("VERGROESZERUNG_Prozeszkraefte_" +
                    Versuchsnummer + "_zoom.png", dpi=300)
        #
        plt.show()
        #
        plt.close()

# Auswertung Ruhekräfte vor der Zerspanung (ab Versuch V0102)
    # HK, Sa, 06.06.2020
    # Vorschub = float(f.split()[0]) # [mm]
    # vSchnitt = float(vc.split()[0])/60.0 # [m/s]
    # D_Zyl = 0.060 # [m] Durchschnittsdurchmesser zur Abschätzung der Zeit bis zur eigentlichen Zerspanung 1mm Weg mit f bis erster Eingriff (abgerundet, konservativ)
    # n= 1.0/Vorschub # [mm/mm] Umdrehungen bei aktuellem Vorschub für 1mm Weg

    # Startzeit für Leerlaufkraftermittlung
    tStart = float(Auswertegrenzen[int(Versuchsnummer) - 1, 6])
    # Endzeit für Leerlaufkraftermittlung
    tEnde = float(Auswertegrenzen[int(Versuchsnummer) - 1, 7])

    Idx_tStart = 0  # Startindex der Auswertung der Leerlaufkräfte
    # if (int(Versuchsnummer) == 1): Idx_tStart = np.argmax(Daten[:,0] > 0.25)
    # # im ersten Versuch ab 0,25s Leerlaufkraft auswerten

    # Endindex der Auswertung der Leerlaufkräfte
    #tmax = np.pi * D_Zyl / vSchnitt * n
    # if (int(Versuchsnummer) < 102): # Normalerweise eine halbe Sekunde
    # Wartezeit vor dem eigentlichen Schnitt (zumindest zum Ende der Versuche
    # ohne Studenten)

    tmax = 12

    if (tStart >= 0 and tEnde >= 0):
        Idx_tStart = np.argmax(Daten[:, 0] > tStart)
        Idx_tmax = np.argmax(Daten[:, 0] > tEnde)
        if (Idx_tmax == 0.0):
            Idx_tmax = len(Daten[:, 0]) - 1
    else:
        Idx_tmax = np.argmax(Daten[:, 0] > tmax)

    FfLeerAvg[int(Versuchsnummer) - 1] = np.nanmean(Ff[Idx_tStart:Idx_tmax])
    FcLeerAvg[int(Versuchsnummer) - 1] = np.nanmean(Fc[Idx_tStart:Idx_tmax])
    FfLeerStd[int(Versuchsnummer) - 1] = np.nanstd(Ff[Idx_tStart:Idx_tmax])
    FcLeerStd[int(Versuchsnummer) - 1] = np.nanstd(Fc[Idx_tStart:Idx_tmax])

    if (DiagrammLeerlaufkraft != 0):
        plt.plot(Daten[Idx_tStart:Idx_tmax, 0],
                 Fc[Idx_tStart:Idx_tmax], label=r'Cutting force')
        plt.plot(Daten[Idx_tStart:Idx_tmax, 0],
                 Ff[Idx_tStart:Idx_tmax], label=r'Feed force')
        LinieFcLeerAvg = [
            FfLeerAvg[int(Versuchsnummer) - 1], FfLeerAvg[int(Versuchsnummer) - 1]]
        LinieFfLeerAvg = [
            FcLeerAvg[int(Versuchsnummer) - 1], FcLeerAvg[int(Versuchsnummer) - 1]]
        plt.plot([Daten[Idx_tStart, 0], Daten[Idx_tmax, 0]], LinieFcLeerAvg, color="darkblue", linewidth=0.5, linestyle="--", label=r'Idle force Fc\_average: ' +
                 str(round(FcLeerAvg[int(Versuchsnummer) - 1], 1)) + '[N], Sigma: ' + str(round(FcLeerStd[int(Versuchsnummer) - 1], 1)))
        plt.plot([Daten[Idx_tStart, 0], Daten[Idx_tmax, 0]], LinieFfLeerAvg, color="darkblue", linewidth=0.5, linestyle="--", label=r'Idle force Ff\_average: ' +
                 str(round(FfLeerAvg[int(Versuchsnummer) - 1], 1)) + '[N], Sigma: ' + str(round(FfLeerStd[int(Versuchsnummer) - 1], 1)))
        plt.xlabel(r'Time [$s$]')
        plt.ylabel(r'Force [N]')
        #plt.title(r'Leerlaufkräfte , Experiment ' + str(Versuchsnummer) + ', vc= ' + vc + ', f= ' + f)
        plt.title(r'Idle forces , Experiment ' + str(Versuchsnummer))
        plt.legend()

        plt.savefig("LEERLAUFKRAEFTE_Prozeszkraefte_" +
                    Versuchsnummer + "_Leerlauf_vor_Versuch.png", dpi=300)
        plt.show()
        plt.close()

    #
#np.savetxt("Prozeszkraefte_gemittelt.txt", (np.transpose(Versuche), np.transpose(Schnittkraft), np.transpose(Passivkraft)), fmt="%2.3f", delimiter=",")
#np.savetxt("Prozeszkraefte_gemittelt.txt", (Versuche, Schnittkraft, Passivkraft), fmt="%2.3f", delimiter=",")
#
#
# plt.plot(Fc_Leerlauf_vorne[0:100])
# plt.plot(Fc_Leerlauf_vorne)
# plt.plot(Ff_Leerlauf_vorne)
# plt.plot(Fc_Leerlauf_hinten)
# plt.plot(Ff_Leerlauf_hinten)
#plt.ylim(-3, 1)

AvgFcLeerAlle = np.nanmean(FcLeerAvg)
AvgFfLeerAlle = np.nanmean(FfLeerAvg)
#
plt.plot(FcLeerAvg, label=r'Cutting force empty')
plt.plot(FfLeerAvg, label=r'Feed force empty')
plt.xlim(0, len(FcLeerAvg) - 1)
plt.xlabel(r'Experiment NO')
plt.ylabel(r'Force [N]')
plt.title(r'Idle forces before experiment; Fc AvgAlle:' +
          str(round(AvgFcLeerAlle, 1)) + 'Ff AvgAlle:' + str(round(AvgFfLeerAlle, 1)))
plt.legend()
#plt.ylim(-3, 1)

plt.savefig(
    "LEERLAUFKRAEFTE_Leerlaufprozeszkraefte_durchschnittlich_vor_Versuch.png", dpi=300)
plt.show()
plt.close()

for i in range(len(FfLeerAvg)):
    if (i > 0):
        if np.isnan(FfLeerAvg[i]):
            FfLeerAvg[i] = FfLeerAvg[i - 1]
        if np.isnan(FcLeerAvg[i]):
            FcLeerAvg[i] = FcLeerAvg[i - 1]

AvgFcLeerAlle = np.nanmean(FcLeerAvg)
AvgFfLeerAlle = np.nanmean(FfLeerAvg)
#
plt.plot(FcLeerAvg, label=r'Schnittkraft leer')
plt.plot(FfLeerAvg, label=r'Vorschubkraft leer')
plt.xlim(0, len(FcLeerAvg) - 1)
plt.xlabel(r'Versuchsnummer')
plt.ylabel(r'Kraft [N]')
plt.title(r'Leerlaufkräfte, vor eigentlichem Versuch (NAN entfernt); Fc AvgAlle:' +
          str(round(AvgFcLeerAlle, 1)) + 'Ff AvgAlle:' + str(round(AvgFfLeerAlle, 1)))
plt.legend()
#plt.ylim(-3, 1)

plt.savefig(
    "LEERLAUFKRAEFTE_Leerlaufprozeszkraefte_durchschnittlich_vor_Versuch_NANentfernt.png", dpi=300)
plt.show()
plt.close()

ProzKraftAusgabe = "Prozeszkraefte_gemittelt.txt"

Datei = open(ProzKraftAusgabe, 'w')
Datei.write('Gemittelte Prozeßkräfte\n')
Datei.write('Die Leerlaufkräfte sind schon auf 1mm Schnittbreite bezogen, bedürfen also keiner Anpassung. Die Standardabweichungen der Prozeßkräfte sind im Prinzip von der Leerlaufkraft unabhängig.\n')
Datei.write('Versuch, Schnittkraft F_c + Leerlauf F_c [N/mm], StdAbw F_c + Leerlauf F_c [N/mm], Vorschubkraft F_f + Leerlauf F_c [N/mm], StdAbw F_f + Leerlauf F_f [N/mm], Fc_LeerAvg [N/mm], Fc_LeerStd [N/mm], Ff_LeerAvg [N/mm], Ff_LeerStd [N/mm], Schnittkraft F_c [N/mm], Vorschubkraft F_f [N/mm\n')
#
for i in range(len(Versuche)):
    # Datei.write(str(Versuche[i]) + ', ' + str(Schnittkraft[i]) + ', ' +
    # str(SigmaSchnittkraft[i])+ ', ' + str(Vorschubkraft[i]) + ', ' +
    # str(SigmaVorschubkraft[i]) + ', ' + str(FcLeerAvg[i]/Wandstaerke) + ', '
    # + str(FcLeerStd[i]/Wandstaerke) + ', ' + str(FfLeerAvg[i]/Wandstaerke) +
    # ', ' + str(FfLeerStd[i]/Wandstaerke) + ', ' +
    # str(Schnittkraft[i]-FcLeerAvg[i]/Wandstaerke) + ', ' +
    # str(Vorschubkraft[i]-FfLeerAvg[i]/Wandstaerke)) # ab Fr, 05.06.2020
    Datei.write(str(Versuche[i]) + ', ' + str(Schnittkraft[i]) + ', ' + str(SigmaSchnittkraft[i]) + ', ' + str(Vorschubkraft[i]) + ', ' + str(SigmaVorschubkraft[i]) + ', ' + str(FcLeerAvg[i] / Wandstaerke) + ', ' + str(FcLeerStd[i] / Wandstaerke) + ', ' + str(
        FfLeerAvg[i] / Wandstaerke) + ', ' + str(FfLeerStd[i] / Wandstaerke) + ', ' + str((Schnittkraft[i] != 0) * (Schnittkraft[i] - FcLeerAvg[i] / Wandstaerke)) + ', ' + str((Vorschubkraft[i] != 0) * (Vorschubkraft[i] - FfLeerAvg[i] / Wandstaerke)))  # ab Fr, 05.06.2020
    Datei.write('\n')
    #~np.isnan(np.nan)*5 + np.nan
#
Datei.close()
