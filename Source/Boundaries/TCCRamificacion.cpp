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
#include <vector>
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

	FGInit = false;
	FGRhoY = NULL;
	FGNx = NULL;
	FGNy = NULL;
	FGRho = FGMx = FGMy = FGE = FGVol = 0.;
	FTiempoActual = 0.;
	FTiempoAnterior = 0.;   // first CalculaCondicionContorno takes DeltaT = Time - FTiempoAnterior
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

// Ghost Junction Method (Hong & Kim 2011), conservative form.
// The ghost junction cell is a finite volume whose faces are the branches' boundary faces. Each branch pipe
// computes the flux through that face itself (RoeM between its end node, which holds the ghost state written
// here, and its first interior cell, every RK stage) and accumulates it; the ghost cell is advanced with exactly
// those amounts (TTubo::TakeEndFlux), so mass, momentum, energy and species are exchanged with one flux per face
// ("the interface flux is used as the boundary flux of neighbouring branches", H&K 2.1.1). Its volume is the
// half-cells the branches give up at their ghost ends. After the update, every branch end node receives the
// ghost state with the normal velocity reconstructed by the scaling function G (H&K Eq.35/37), using the velocity
// of the branch's adjacent interior cell. Equally-spaced branch normals (OpenWAM has no branch angles).
void TCCRamificacion::CalculaCondicionContornoGJM(double DeltaT) {
	const int N = FNumeroTubosCC;
	const int nsp = FNumeroEspecies - FIntEGR;
	const double ARef = __cons::ARef;
	const double TWOPI = 6.283185307179586;
	const int MAXB = 16;
	if(N > MAXB)
		throw Exception("GJM: too many branches at junction");

	int sgn[MAXB], endp[MAXB];   // sgn +1: left pipe end (N_i = pipe +x), -1: right end (N_i = pipe -x)
	double Ubn[MAXB], acell[MAXB], gam[MAXB], Acara[MAXB];
	for(int i = 0; i < N; i++) {
		TTubo* P = FTuboExtremo[i].Pipe;
		sgn[i] = (FIndiceCC[i] == 0) ? 1 : -1;
		endp[i] = (FIndiceCC[i] == 0) ? 0 : 1;
		const int c = (endp[i] == 0) ? 1 : P->getNin() - 2;   // adjacent interior cell
		Ubn[i] = sgn[i] * P->GetVelocidad(c) * ARef;          // its velocity along N_i (junction -> pipe)
		acell[i] = P->GetAsonido(c) * ARef;
		gam[i] = P->GetGamma(FNodoFin[i]);
		Acara[i] = P->GetAreaCara(endp[i]) * P->getNumeroConductos();
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
	std::vector<double> Y(nsp);
	if(!FGInit) {                                             // initial ghost: area-weighted end states, at rest
		double Aw = 0., rw = 0., pw = 0., gw = 0., vsum = 0.;
		for(int i = 0; i < N; i++) {
			TTubo* P = FTuboExtremo[i].Pipe;
			int nd = FNodoFin[i];
			double A = P->GetArea(nd) * P->getNumeroConductos();
			vsum += 0.5 * A * P->getXRef();                    // the half-cell this branch end gives up
			Aw += A;
			rw += P->GetDensidad(nd) * A;
			pw += __units::BarToPa(P->GetPresion(nd)) * A;
			gw += gam[i] * A;
		}
		FGVol = vsum;
		FGRho = rw / Aw;
		FGMx = 0.;
		FGMy = 0.;
		FGE = (pw / Aw) / (gw / Aw - 1.0);
		for(int k = 0; k < nsp; k++) {
			double yw = 0.;
			for(int i = 0; i < N; i++)
				yw += FTuboExtremo[i].Pipe->GetFraccionMasicaCC(FIndiceCC[i], k) * FTuboExtremo[i].Pipe->GetDensidad(FNodoFin[i]) *
					  FTuboExtremo[i].Pipe->GetArea(FNodoFin[i]) * FTuboExtremo[i].Pipe->getNumeroConductos();
			FGRhoY[k] = yw / Aw;
		}
		for(int i = 0; i < N; i++) {                          // anything already exchanged belongs to the ghost cell
			double m, q, e;
			FTuboExtremo[i].Pipe->TakeEndFlux(endp[i], m, q, e, &Y[0]);
			FGRho -= sgn[i] * m / FGVol;
			FGE -= sgn[i] * e / FGVol;
			for(int k = 0; k < nsp; k++)
				FGRhoY[k] -= sgn[i] * Y[k] / FGVol;
		}
		FGInit = true;
	} else {
		// Ghost state at the start of the step (the state the branches used for their face fluxes).
		double gg = 0.;
		for(int i = 0; i < N; i++)
			gg += gam[i];
		gg /= N;
		const double ug = FGMx / FGRho, vg = FGMy / FGRho;
		double pg = (gg - 1.0) * (FGE - 0.5 * FGRho * (ug * ug + vg * vg));
		if(pg < 1.0)
			pg = 1.0;
		// Advance the ghost cell with exactly what crossed each branch face this step.
		double dM = 0., dMx = 0., dMy = 0., dE = 0.;
		std::vector<double> dY(nsp, 0.);
		for(int i = 0; i < N; i++) {
			double m, q, e;
			FTuboExtremo[i].Pipe->TakeEndFlux(endp[i], m, q, e, &Y[0]);
			const double mout = sgn[i] * m;   // mass leaving the ghost through face i (along N_i)
			const double Nx = FGNx[i], Ny = FGNy[i];
			const double Un_g = ug * Nx + vg * Ny;
			// Normal-momentum impulse through the face: (rho u^2 + p) A dt, identical in the N_i frame for both
			// end orientations; tangential ghost momentum is carried out only by outflow (H&K Eq.17).
			const double Utx = (mout > 0.) ? ug - Un_g * Nx : 0.;
			const double Uty = (mout > 0.) ? vg - Un_g * Ny : 0.;
			dM += mout;
			dE += sgn[i] * e;
			dMx += q * Nx + mout * Utx - pg * Nx * Acara[i] * DeltaT;   // minus the wall reaction (H&K Eq.11)
			dMy += q * Ny + mout * Uty - pg * Ny * Acara[i] * DeltaT;
			for(int k = 0; k < nsp; k++)
				dY[k] += sgn[i] * Y[k];
		}
		FGRho -= dM / FGVol;
		FGMx -= dMx / FGVol;
		FGMy -= dMy / FGVol;
		FGE -= dE / FGVol;
		for(int k = 0; k < nsp; k++)
			FGRhoY[k] -= dY[k] / FGVol;
		if(FGRho < 1e-6) {
			printf("WARNING: GJM junction %d: ghost density %g clamped (junction volume too small for the time step?)\n", FNumeroCC,
				   FGRho);
			FGRho = 1e-6;
		}
	}
	for(int k = 0; k < nsp; k++) {                            // junction outflow composition
		FFraccionMasicaEspecie[k] = FGRhoY[k] / FGRho;
		if(FFraccionMasicaEspecie[k] < 0.)
			FFraccionMasicaEspecie[k] = 0.;
	}

	// New ghost state, written as each branch's end-node (ghost) state for the next step's face flux.
	double gg = 0.;
	for(int i = 0; i < N; i++)
		gg += gam[i];
	gg /= N;
	const double ug = FGMx / FGRho, vg = FGMy / FGRho;
	double pg = (gg - 1.0) * (FGE - 0.5 * FGRho * (ug * ug + vg * vg));
	if(pg < 1.0)
		pg = 1.0;
	for(int i = 0; i < N; i++) {
		const double Un_g = ug * FGNx[i] + vg * FGNy[i];
		const double dUn = Ubn[i] - Un_g;
		const double delta = (Ubn[i] < 0.0) ? 0.0 : 1.0;      // 0 inflow to junction, 1 outflow (Eq 35)
		const double G = delta * fmin(fabs(dUn) / acell[i], 1.0);
		const double un = (Ubn[i] - G * dUn) * sgn[i];        // reconstructed ghost normal velocity, pipe +x (Eq 37)
		// Use the branch's own gamma so the pipe recovers exactly rho = FGRho, p = pg from Landa/Beta/Entropia.
		const double an = sqrt(gam[i] * pg / FGRho);
		const double adim = an / ARef, vdim = un / ARef, pbar = __units::PaToBar(pg);
		const double G3 = __Gamma::G3(gam[i]), G5 = __Gamma::G5(gam[i]);
		FTuboExtremo[i].Landa = adim + G3 * vdim;
		FTuboExtremo[i].Beta = adim - G3 * vdim;
		FTuboExtremo[i].Entropia = adim / pow(pbar, G5);
	}
}

double TCCRamificacion::getGhostMass() const {
	return FGInit ? FGRho * FGVol : 0.;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TCCRamificacion::CalculaCondicionContorno(double Time) {
	try {
		double DeltaT;

		FTiempoActual = Time;
		DeltaT = FTiempoActual - FTiempoAnterior;
		FTiempoAnterior = FTiempoActual;

		// GJM (Ghost Junction Method): the production junction, run once per step (first call, DeltaT>0),
		// updating every branch's Landa/Beta/Entropia + species; later per-pipe calls this step return.
		if(DeltaT > 1e-12)
			CalculaCondicionContornoGJM(DeltaT);

	} catch(exception & N) {
		std::cout << "ERROR: TCCRamificacion::CalculaCondicionContorno en la condicion de contorno: " << FNumeroCC << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#pragma package(smart_init)
