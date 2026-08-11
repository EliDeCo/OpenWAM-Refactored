/* --------------------------------------------------------------------------------*\
==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 \\ |  X  | //  W ave     |
 \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*-------------------------------------------------------------------------------- */

// ---------------------------------------------------------------------------
#pragma hdrstop

#include "TCCRamificacion.h"
//#include <cmath>
#include <iostream>
#include "TTubo.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

TCCRamificacion::TCCRamificacion(nmTypeBC TipoCC, int numCC, nmTipoCalculoEspecies SpeciesModel, int numeroespecies,
								 nmCalculoGamma GammaCalculation, bool ThereIsEGR) :
	TCondicionContorno(TipoCC, numCC, SpeciesModel, numeroespecies, GammaCalculation, ThereIsEGR) {

	FTuboExtremo = NULL;

	FNodoFin = NULL;
	FIndiceCC = NULL;
	FEntropia = NULL;
	FSeccionTubo = NULL;
	FVelocity = NULL;
	FDensidad = NULL;
	FNumeroTubo = NULL;

	FCC = NULL;
	FCD = NULL;

	FMasaEspecie = NULL;

	FUseGJM = -1;
	FGInit = false;
	FGRhoY = NULL;
	FGNx = NULL;
	FGNy = NULL;
	FGRho = FGMx = FGMy = FGE = FGVol = 0.;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

TCCRamificacion::~TCCRamificacion() {

	if(FTuboExtremo != NULL)
		delete[] FTuboExtremo;
	if(FNodoFin != NULL)
		delete[] FNodoFin;
	if(FIndiceCC != NULL)
		delete[] FIndiceCC;
	if(FEntropia != NULL)
		delete[] FEntropia;
	if(FSeccionTubo != NULL)
		delete[] FSeccionTubo;
	if(FVelocity != NULL)
		delete[] FVelocity;
	if(FDensidad != NULL)
		delete[] FDensidad;
	if(FNumeroTubo != NULL)
		delete[] FNumeroTubo;

	if(FCC != NULL)
		delete[] FCC;
	if(FCD != NULL)
		delete[] FCD;

	if(FMasaEspecie != NULL)
		delete[] FMasaEspecie;

	if(FGRhoY != NULL)
		delete[] FGRhoY;
	if(FGNx != NULL)
		delete[] FGNx;
	if(FGNy != NULL)
		delete[] FGNy;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TCCRamificacion::AsignaTubos(int NumberOfPipes, TTubo **Pipe) {
	try {
		int i = 0;
		int ContadorTubosRamificacion = 0;

		ContadorTubosRamificacion = 0;

		for(int i = 0; i < NumberOfPipes; i++) {
			if(Pipe[i]->getNodoIzq() == FNumeroCC || Pipe[i]->getNodoDer() == FNumeroCC) {
				ContadorTubosRamificacion++;
			}
		}

		FTuboExtremo = new stTuboExtremo[ContadorTubosRamificacion];
		FNodoFin = new int[ContadorTubosRamificacion];
		FIndiceCC = new int[ContadorTubosRamificacion];
		FCC = new double*[ContadorTubosRamificacion];
		FCD = new double*[ContadorTubosRamificacion];
		FEntropia = new double[ContadorTubosRamificacion];
		FSeccionTubo = new double[ContadorTubosRamificacion];
		FVelocity = new double[ContadorTubosRamificacion];
		FDensidad = new double[ContadorTubosRamificacion];
		FNumeroTubo = new int[ContadorTubosRamificacion];

		for(int i = 0; i < ContadorTubosRamificacion; i++) {
			FTuboExtremo[i].Pipe = NULL;
			FVelocity[i] = 0;
		}

		while(FNumeroTubosCC < ContadorTubosRamificacion && i < NumberOfPipes) {
			if(Pipe[i]->getNodoIzq() == FNumeroCC) {
				FTuboExtremo[FNumeroTubosCC].Pipe = Pipe[i];
				FTuboExtremo[FNumeroTubosCC].TipoExtremo = nmLeft;
				FNodoFin[FNumeroTubosCC] = 0;
				FIndiceCC[FNumeroTubosCC] = 0;
				FNumeroTubo[FNumeroTubosCC] = Pipe[i]->getNumeroTubo() - 1;
				FCC[FNumeroTubosCC] = &(FTuboExtremo[FNumeroTubosCC].Beta);
				FCD[FNumeroTubosCC] = &(FTuboExtremo[FNumeroTubosCC].Landa);
				FSeccionTubo[FNumeroTubosCC] = __geom::Circle_area(Pipe[i]->GetDiametro(FNodoFin[FNumeroTubosCC]));
				FNumeroTubosCC++;
			}
			if(Pipe[i]->getNodoDer() == FNumeroCC) {
				FTuboExtremo[FNumeroTubosCC].Pipe = Pipe[i];
				FTuboExtremo[FNumeroTubosCC].TipoExtremo = nmRight;
				FNodoFin[FNumeroTubosCC] = Pipe[i]->getNin() - 1;
				FIndiceCC[FNumeroTubosCC] = 1;
				FNumeroTubo[FNumeroTubosCC] = Pipe[i]->getNumeroTubo() - 1;
				FCC[FNumeroTubosCC] = &(FTuboExtremo[FNumeroTubosCC].Landa);
				FCD[FNumeroTubosCC] = &(FTuboExtremo[FNumeroTubosCC].Beta);
				FSeccionTubo[FNumeroTubosCC] = __geom::Circle_area(Pipe[i]->GetDiametro(FNodoFin[FNumeroTubosCC]));
				FNumeroTubosCC++;
			}
			i++;
		}

		// Inicializacion del transporte de especies quimicas.
		FFraccionMasicaEspecie = new double[FNumeroEspecies - FIntEGR];
		FMasaEspecie = new double[FNumeroEspecies - FIntEGR];
		for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
			FFraccionMasicaEspecie[i] = FTuboExtremo[0].Pipe->GetFraccionMasicaInicial(i);
			// Se inicializa con el Pipe 0 de modo arbitrario.
		}

	} catch(exception & N) {
		std::cout << "ERROR: TCCRamificacion::AsignaTubos en la condicion de contorno: " << FNumeroCC << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TCCRamificacion::TuboCalculandose(int TuboActual) {
	try {
		FTuboActual = TuboActual;
		if(FTuboActual == 10000) {
			FTiempoActual = FTuboExtremo[0].Pipe->getTime1();
		} else {
			for(int i = 0; i < FNumeroTubosCC; i++) {
				if(FNumeroTubo[i] == FTuboActual) {
					FTiempoActual = FTuboExtremo[i].Pipe->getTime1();
				}
			}
		}
	} catch(exception & N) {
		std::cout << "ERROR: TCCRamificacion::TuboCalculandose en la condicion de contorno: " << FNumeroCC << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// GJM: quasi-1D RoeM interface flux (Hong & Kim 2011, Eq 40), SI units. Returns the physical flux
// per unit area: frho, frhoun (= rho*U^2+p momentum), frhoet. Mirrors TTubo::RoeMFlux.
void TCCRamificacion::RoeM1DFlux(double rL, double uL, double pL, double gL, double rR, double uR, double pR,
								 double gR, double& frho, double& frhoun, double& frhoet) {
	const double G1L = gL - 1.0, G1R = gR - 1.0;
	const double HL = gL * pL / (G1L * rL) + 0.5 * uL * uL;
	const double HR = gR * pR / (G1R * rR) + 0.5 * uR * uR;
	const double FL0 = rL * uL, FL1 = rL * uL * uL + pL, FL2 = rL * uL * HL;
	const double FR0 = rR * uR, FR1 = rR * uR * uR + pR, FR2 = rR * uR * HR;
	const double srL = sqrt(rL), srR = sqrt(rR), sden = srL + srR;
	const double rhat = srL * srR;
	const double uhat = (srL * uL + srR * uR) / sden;
	const double Hhat = (srL * HL + srR * HR) / sden;
	const double G1 = (srL * G1L + srR * G1R) / sden;
	const double ahat = sqrt(G1 * (Hhat - 0.5 * uhat * uhat));
	const double Mhat = uhat / ahat;
	const double b1 = fmax(0.0, fmax(uhat + ahat, uR + ahat));
	const double b2 = fmin(0.0, fmin(uhat - ahat, uL - ahat));
	const double prat = fmin(pL / pR, pR / pL);
	const double hexp = 1.0 - prat;
	const double fc = (fabs(uhat) < 1.0e-30) ? 1.0 : pow(fabs(Mhat), hexp);
	const double dQ0 = rR - rL, dQ1 = rR * uR - rL * uL, dQ2 = rR * HR - rL * HL;
	const double dp = pR - pL, dH = HR - HL;
	const double sfac = (rR - rL) - fc * dp / (ahat * ahat);
	const double B0 = sfac, B1 = sfac * uhat, B2 = sfac * Hhat + rhat * dH;
	const double bden = b1 - b2, w = b1 * b2 / bden, wB = w / (1.0 + fabs(Mhat));
	frho = (b1 * FL0 - b2 * FR0) / bden + w * dQ0 - wB * B0;
	frhoun = (b1 * FL1 - b2 * FR1) / bden + w * dQ1 - wB * B1;
	frhoet = (b1 * FL2 - b2 * FR2) / bden + w * dQ2 - wB * B2;
}

// Ghost Junction Method (Hong & Kim 2011): persistent 2-D ghost cell (SI), RoeM interface flux +
// scaling function G (Eq 35/37), equally-spaced branch normals (user ruling: OpenWAM has no angles).
// Runs once per step and writes each branch's new (Landa,Beta,Entropia). See Notes/GJM-Implementation-Recipe.md.
void TCCRamificacion::CalculaCondicionContornoGJM(double DeltaT) {
	const int N = FNumeroTubosCC;
	const int nsp = FNumeroEspecies - FIntEGR;
	const double ARef = __cons::ARef;
	const double TWOPI = 6.283185307179586;
	const int MAXB = 16;
	if(N > MAXB)
		throw Exception("GJM: too many branches at junction");

	double rho[MAXB], up[MAXB], Ubn[MAXB], pr[MAXB], aa[MAXB], gam[MAXB], Ar[MAXB];
	int sgn[MAXB];   // +1 left end (N_i = pipe +x), -1 right end (N_i = pipe -x)
	for(int i = 0; i < N; i++) {
		TTubo* P = FTuboExtremo[i].Pipe;
		int nd = FNodoFin[i];
		rho[i] = P->GetDensidad(nd);
		up[i] = P->GetVelocidad(nd) * ARef;                   // pipe +x velocity (SI)
		pr[i] = __units::BarToPa(P->GetPresion(nd));
		aa[i] = P->GetAsonido(nd) * ARef;
		gam[i] = P->GetGamma(nd);
		Ar[i] = FSeccionTubo[i];
		sgn[i] = (FIndiceCC[i] == 0) ? 1 : -1;
		Ubn[i] = up[i] * sgn[i];                              // velocity along N_i (junction -> pipe)
	}

	if(FGNx == NULL) {                                        // equally-spaced branch normals + species buffer
		FGNx = new double[N];
		FGNy = new double[N];
		FGRhoY = new double[nsp];
		for(int i = 0; i < N; i++) {
			double th = TWOPI * i / N;
			FGNx[i] = cos(th);
			FGNy[i] = sin(th);
		}
	}
	if(!FGInit) {                                             // init ghost = area-weighted, at rest; skip update this step
		double Aw = 0., rw = 0., pw = 0., gw = 0., vsum = 0.;
		for(int i = 0; i < N; i++) {
			TTubo* P = FTuboExtremo[i].Pipe;
			double dx = P->getLongitudTotal() / (P->getNin() - 1);
			vsum += Ar[i] * dx;
			Aw += Ar[i];
			rw += rho[i] * Ar[i];
			pw += pr[i] * Ar[i];
			gw += gam[i] * Ar[i];
		}
		FGVol = vsum / N;
		FGRho = rw / Aw;
		double gg0 = gw / Aw, pg0 = pw / Aw;
		FGMx = 0.;
		FGMy = 0.;
		FGE = pg0 / (gg0 - 1.0);
		for(int k = 0; k < nsp; k++) {
			double yw = 0.;
			for(int i = 0; i < N; i++)
				yw += FTuboExtremo[i].Pipe->GetFraccionMasicaCC(FIndiceCC[i], k) * rho[i] * Ar[i];
			FGRhoY[k] = yw / Aw;
		}
		FGInit = true;
		return;
	}

	const double ug = FGMx / FGRho, vg = FGMy / FGRho;        // ghost primitives
	double gg = 0.;
	for(int i = 0; i < N; i++)
		gg += gam[i];
	gg /= N;
	double pg = (gg - 1.0) * (FGE - 0.5 * FGRho * (ug * ug + vg * vg));
	if(pg < 1.0)
		pg = 1.0;

	double UnLrec[MAXB], sumY[MAXB];
	double sumR = 0., sumMx = 0., sumMy = 0., sumE = 0., sumWx = 0., sumWy = 0.;
	for(int k = 0; k < nsp; k++)
		sumY[k] = 0.;
	for(int i = 0; i < N; i++) {                              // interface fluxes (N_i frame) + ghost accumulation
		const double Nx = FGNx[i], Ny = FGNy[i];
		const double Un_g = ug * Nx + vg * Ny;
		const double Utx = ug - Un_g * Nx, Uty = vg - Un_g * Ny;
		const double dUn = Ubn[i] - Un_g;
		const double delta = (Ubn[i] < 0.0) ? 0.0 : 1.0;      // 0 inflow to junction, 1 outflow (Eq 35)
		const double G = delta * fmin(fabs(dUn) / aa[i], 1.0);
		UnLrec[i] = Ubn[i] - G * dUn;                         // reconstructed ghost normal velocity (Eq 37)

		double fr, fmn, fe;
		RoeM1DFlux(FGRho, UnLrec[i], pg, gg, rho[i], Ubn[i], pr[i], gam[i], fr, fmn, fe);
		const double Utxu = (fr >= 0.0) ? Utx : 0.0;          // tangential advected only from the ghost side
		const double Utyu = (fr >= 0.0) ? Uty : 0.0;
		sumR += fr * Ar[i];
		sumE += fe * Ar[i];
		sumMx += (fmn * Nx + fr * Utxu) * Ar[i];
		sumMy += (fmn * Ny + fr * Utyu) * Ar[i];
		sumWx += pg * Nx * Ar[i];                             // wall reaction (Eq 11)
		sumWy += pg * Ny * Ar[i];
		for(int k = 0; k < nsp; k++) {
			double Yup = (fr >= 0.0) ? (FGRhoY[k] / FGRho) : FTuboExtremo[i].Pipe->GetFraccionMasicaCC(FIndiceCC[i], k);
			sumY[k] += fr * Yup * Ar[i];
		}
	}

	const double iv = DeltaT / FGVol;                         // ghost cell update (Eq 7a + wall force)
	FGRho -= iv * sumR;
	FGMx -= iv * (sumMx - sumWx);
	FGMy -= iv * (sumMy - sumWy);
	FGE -= iv * sumE;
	for(int k = 0; k < nsp; k++)
		FGRhoY[k] -= iv * sumY[k];
	if(FGRho < 1e-6)
		FGRho = 1e-6;
	if(getenv("OPENWAM_GJM_DIAG") != NULL) {
		// netMass = net mass flux out of ghost, netH0 = net stagnation-enthalpy flux out of ghost (Sum mdot*h0);
		// thru = enthalpy-flux throughput. A conservative junction absorbs the imbalance in the ghost cell,
		// which returns to steady (no cumulative insertion). See ramification-junction-energy-nonconservation.
		double thru = 0.;
		for(int i = 0; i < N; i++) {
			double frr, fmm, fee;
			double Un_g = ug * FGNx[i] + vg * FGNy[i];
			RoeM1DFlux(FGRho, Ubn[i] - ((Ubn[i] < 0.) ? 0. : fmin(fabs(Ubn[i] - Un_g) / aa[i], 1.)) * (Ubn[i] - Un_g),
					   pg, gg, rho[i], Ubn[i], pr[i], gam[i], frr, fmm, fee);
			thru += fabs(fee * Ar[i]);
		}
		printf("GJMDIAG cc=%d t=%.6e rho_g=%.5f p_g=%.1f E_g=%.2f netMass=%.3e netH0=%.3e relH0=%.3e\n", FNumeroCC,
			   FTiempoActual, FGRho, pg, FGE, sumR, sumE, (thru > 1e-30) ? fabs(sumE) / thru : 0.);
	}
	for(int k = 0; k < nsp; k++) {                            // junction outflow composition
		FFraccionMasicaEspecie[k] = FGRhoY[k] / FGRho;
		if(FFraccionMasicaEspecie[k] < 0.)
			FFraccionMasicaEspecie[k] = 0.;
	}

	for(int i = 0; i < N; i++) {                              // branch boundary-cell update (Eq 18) + write-back
		TTubo* P = FTuboExtremo[i].Pipe;
		int nd = FNodoFin[i], nd1 = nd + sgn[i];
		double dx = P->getLongitudTotal() / (P->getNin() - 1);
		double A0 = P->GetArea(nd), Af = 0.5 * (A0 + P->GetArea(nd1));
		double etv = pr[i] / (gam[i] - 1.0) + 0.5 * rho[i] * up[i] * up[i];
		double Q0 = rho[i] * A0, Q1 = rho[i] * up[i] * A0, Q2 = etv * A0;

		double r1 = P->GetDensidad(nd1), u1 = P->GetVelocidad(nd1) * ARef, p1 = __units::BarToPa(P->GetPresion(nd1)),
			   g1 = P->GetGamma(nd1);
		double Fi0, Fi1, Fi2, Fj0, Fj1, Fj2;
		double ugx = UnLrec[i] * sgn[i];                      // ghost velocity in pipe +x
		if(sgn[i] == 1) {                                     // left end: junction face on the left
			RoeM1DFlux(rho[i], up[i], pr[i], gam[i], r1, u1, p1, g1, Fi0, Fi1, Fi2);
			RoeM1DFlux(FGRho, ugx, pg, gg, rho[i], up[i], pr[i], gam[i], Fj0, Fj1, Fj2);
		} else {                                              // right end: junction face on the right
			RoeM1DFlux(r1, u1, p1, g1, rho[i], up[i], pr[i], gam[i], Fi0, Fi1, Fi2);
			RoeM1DFlux(rho[i], up[i], pr[i], gam[i], FGRho, ugx, pg, gg, Fj0, Fj1, Fj2);
		}
		double c = DeltaT / dx, nQ0, nQ1, nQ2;
		if(sgn[i] == 1) {
			nQ0 = Q0 - c * (Af * Fi0 - A0 * Fj0);
			nQ1 = Q1 - c * (Af * Fi1 - A0 * Fj1);
			nQ2 = Q2 - c * (Af * Fi2 - A0 * Fj2);
		} else {
			nQ0 = Q0 - c * (A0 * Fj0 - Af * Fi0);
			nQ1 = Q1 - c * (A0 * Fj1 - Af * Fi1);
			nQ2 = Q2 - c * (A0 * Fj2 - Af * Fi2);
		}
		double rn = nQ0 / A0;
		if(rn < 1e-6)
			rn = 1e-6;
		double un = nQ1 / nQ0, etn = nQ2 / A0;
		double pn = (gam[i] - 1.0) * (etn - 0.5 * rn * un * un);
		if(pn < 1.0)
			pn = 1.0;
		double an = sqrt(gam[i] * pn / rn);
		double adim = an / ARef, vdim = un / ARef, pbar = __units::PaToBar(pn);
		double G3 = __Gamma::G3(gam[i]), G5 = __Gamma::G5(gam[i]);
		FTuboExtremo[i].Landa = adim + G3 * vdim;
		FTuboExtremo[i].Beta = adim - G3 * vdim;
		FTuboExtremo[i].Entropia = adim / pow(pbar, G5);
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TCCRamificacion::CalculaCondicionContorno(double Time) {
	try {
		double sonido_supuesta_ad, sonido_ant_ad, entropia_entrante, corr_entropia;
		double suma1 = 0., suma2 = 0., sm1 = 0., sm2 = 0., sm3 = 0.;
		int TuboCalculado = 0;
		double DeltaT, MasaTotal = 0., g, m, FraccionMasicaAcum = 0.;
		// Necesarias para el calculo de especies en la BC.

		sonido_supuesta_ad = 0.;

		FTiempoActual = Time;
		DeltaT = FTiempoActual - FTiempoAnterior;
		FTiempoAnterior = FTiempoActual;

		// GJM (Ghost Junction Method) path. Runs once per step (first call, DeltaT>0) and updates ALL
		// branch characteristics; subsequent per-pipe calls this step just return. See CalculaCondicionContornoGJM.
		if(FUseGJM == -1)
			FUseGJM = (getenv("OPENWAM_GJM") != NULL) ? 1 : 0;
		if(FUseGJM == 1) {
			if(DeltaT > 1e-12)
				CalculaCondicionContornoGJM(DeltaT);
			return;
		}

		if(FTuboActual == 10000) {
			TuboCalculado = FTuboActual;
			FGamma = FTuboExtremo[0].Pipe->GetGamma(FNodoFin[0]);
		} else {
			for(int i = 0; i < FNumeroTubosCC; i++) {
				if(FNumeroTubo[i] == FTuboActual) {
					TuboCalculado = i;
				}
			}
			FGamma = FTuboExtremo[TuboCalculado].Pipe->GetGamma(FNodoFin[TuboCalculado]);
		}
		// FGamma=FTuboExtremo[TuboCalculado].Pipe->GetGamma(FNodoFin[TuboCalculado]);
		FGamma1 = __Gamma::G1(FGamma);
		FGamma3 = __Gamma::G3(FGamma);
		FGamma4 = __Gamma::G4(FGamma);

		for(int i = 0; i < FNumeroTubosCC; i++) {
			FEntropia[i] = FTuboExtremo[i].Entropia;
		}

		do {
			/* Determinacion de la velocidad del sonido en la ramificacion. */
			suma1 = 0.;
			suma2 = 0.;
			for(int i = 0; i < FNumeroTubosCC; i++) {
				suma1 = suma1 + (*FCC[i]) * FSeccionTubo[i] / pow2(FEntropia[i]);
				suma2 = suma2 + FTuboExtremo[i].Entropia * FSeccionTubo[i] / pow2(FEntropia[i]);
			}
			sonido_ant_ad = sonido_supuesta_ad;
			sonido_supuesta_ad = suma1 / suma2; // Velocity del sonido adimensionalizada (si las variables fuesen dimensionales).
			// Es una especie de promedio respecto de la entropia de cada tubo.

			/* Determinacion de la velocidad de cada tubo de la ramificacion. Esta
			 velocidad sera positiva si el flujo sale del tubo (tubo saliente) y
			 negativa si el flujo entra en el tubo (tubo entrante). En realidad se trata
			 de la velocidad solo para extremo derecho. Para el extremo izquierdo,esta multiplicada
			 por un signo negativo. */
			for(int i = 0; i < FNumeroTubosCC; i++) {
				FVelocity[i] = (*FCC[i] - sonido_supuesta_ad * FTuboExtremo[i].Entropia) / FGamma3;
			}

			/* Calculo de la entropia de los tubos entrantes (el flujo entra en ellos).
			 Esta entropia es igual para todos ellos y se calcula como un balance del
			 flujo que llega de los tubos salientes (de velocidad positiva). */
			sm1 = 0.;
			sm2 = 0.;
			sm3 = 0.;
			for(int i = 0; i < FNumeroTubosCC; i++) {
				sm3 = sm3 + FTuboExtremo[i].Entropia;
				if(FVelocity[i] > 2e-6) {
					sm1 = sm1 + FVelocity[i] * FSeccionTubo[i] * FEntropia[i];
					sm2 = sm2 + FVelocity[i] * FSeccionTubo[i];
				}
			}

			if(sm2 < 2e-6) {
				entropia_entrante = sm3 / FNumeroTubosCC;
			} else {
				/* Desde el punto de vista teorico esta es la forma correcta. La
				 formula anterior es para evitar errores de indeterminacion si
				 sm2 es muy pequena,por lo que se acepta como aproximacion. */
				entropia_entrante = sm1 / sm2;
			}
			for(int i = 0; i < FNumeroTubosCC; i++) {
				FEntropia[i] = FTuboExtremo[i].Entropia;
				if(FVelocity[i] < 0) {   // Flujo entrante al tubo
					FEntropia[i] = entropia_entrante;
				}
			}
		} while((sonido_supuesta_ad - sonido_ant_ad) / (sonido_ant_ad + 0.01) > 1e-4);

		/* Calculo de las caracteristicas y la entropia en los extremos del tubo que se
		 esta calculando, una vez resuelta la condicion de contorno */
		if(TuboCalculado != 10000) {
			corr_entropia = FTuboExtremo[TuboCalculado].Entropia / FEntropia[TuboCalculado];
			*FCC[TuboCalculado] = (*FCC[TuboCalculado] + FGamma3 * FVelocity[TuboCalculado] * (corr_entropia - 1)) / corr_entropia;
			*FCD[TuboCalculado] = *FCC[TuboCalculado] - FGamma1 * FVelocity[TuboCalculado];
			FTuboExtremo[TuboCalculado].Entropia = FEntropia[TuboCalculado];

			double ason = (*FCC[TuboCalculado] + *FCD[TuboCalculado]) / 2;
			double Machx = fabs(FVelocity[TuboCalculado]) / ason;
			if(Machx > 1) {
				printf("Sonic condition in boundary: %d\n", FNumeroCC);
				// double Machy = Machx / fabs(Machx) * sqrt
				// ((pow(Machx, 2) + 2. / FGamma1) / (FGamma4 * pow(Machx, 2) - 1.));
				// double asonido = (*FCC[TuboCalculado] + *FCD[TuboCalculado]) / 2;
				// double Sonidoy = asonido * sqrt
				// ((FGamma3 * pow(Machx, 2) + 1.) / (FGamma3 * pow(Machy, 2) + 1.));
				//
				// double Velocidady = Sonidoy * Machy;
				ReduceSubsonicFlow(ason, FVelocity[TuboCalculado], FGamma);
				*FCC[TuboCalculado] = ason + FGamma3 * FVelocity[TuboCalculado];
				*FCD[TuboCalculado] = ason - FGamma3 * FVelocity[TuboCalculado];
			}
		} else {
			for(int i = 0; i < FNumeroTubosCC; i++) {
				corr_entropia = FTuboExtremo[i].Entropia / FEntropia[i];
				*FCC[i] = (*FCC[i] + FGamma3 * FVelocity[i] * (corr_entropia - 1)) / corr_entropia;
				*FCD[i] = *FCC[i] - FGamma1 * FVelocity[i];
				FTuboExtremo[i].Entropia = FEntropia[i];
				double Machx = fabs(*FCC[i] - *FCD[i]) / (*FCC[i] + *FCD[i]) * 2 / FGamma1;
				if(Machx > 1) {
					printf("Sonic condition in boundary: %d\n", FNumeroCC);
					double Machy = Machx / fabs(Machx) * sqrt((pow2(Machx) + 2. / FGamma1) / (FGamma4 * pow2(Machx) - 1.));
					double asonido = (*FCC[i] + *FCD[i]) / 2;
					double Sonidoy = asonido * sqrt((FGamma3 * pow2(Machx) + 1.) / (FGamma3 * pow2(Machy) + 1.));

					double Velocidady = Sonidoy * Machy;
					*FCC[i] = Sonidoy + FGamma3 * Velocidady;
					*FCD[i] = Sonidoy - FGamma3 * Velocidady;
				}
			}
		}

		// Transporte de especies quimicas.
		for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
			FMasaEspecie[j] = 0.;
		}
		for(int i = 0; i < FNumeroTubosCC; i++) {
			if(FVelocity[i] > 0.) {   // Flujo Saliente del tubo
				FDensidad[i] = pow(((*FCC[i] + *FCD[i]) / 2) / FTuboExtremo[i].Entropia, FGamma4);
				g = FDensidad[i] * FSeccionTubo[i] * FVelocity[i];
				m = g * DeltaT;
				MasaTotal += m;
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					FMasaEspecie[j] += FTuboExtremo[i].Pipe->GetFraccionMasicaCC(FIndiceCC[i], j) * m;
				}
			}
		}

		if(MasaTotal != 0) {
			for(int j = 0; j < FNumeroEspecies - 2; j++) {
				FFraccionMasicaEspecie[j] = FMasaEspecie[j] / MasaTotal;
				FraccionMasicaAcum += FFraccionMasicaEspecie[j];
			}
			FFraccionMasicaEspecie[FNumeroEspecies - 2] = 1. - FraccionMasicaAcum;
			if(FHayEGR)
				FFraccionMasicaEspecie[FNumeroEspecies - 1] = FMasaEspecie[FNumeroEspecies - 1] / MasaTotal;
		}

	} catch(exception & N) {
		std::cout << "ERROR: TCCRamificacion::CalculaCondicionContorno en la condicion de contorno: " << FNumeroCC << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#pragma package(smart_init)
