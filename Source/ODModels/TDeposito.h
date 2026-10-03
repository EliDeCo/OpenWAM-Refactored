/* --------------------------------------------------------------------------------*\
==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
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
#ifndef TDepositoH
#define TDepositoH
#include <iostream>
#ifdef __BORLANDC__
#include <vcl.h>
#endif
//#include <cmath>
#include <vector>
#include <cstdio>

#include "Globales.h"

class TCondicionContorno;
class TCompresor;

class TDeposito {
  protected:
	// ---------------------------------------------------------------------------
	// VARIABLES PRIVADAS
	// ---------------------------------------------------------------------------

	TCondicionContorno **FCCDeposito;
	TCondicionContorno **FCCUnionEntreDep;

	// Constantes termodinamicas

	double FPresRef;

	// Variables geometricas

	int FNumeroDeposito;
	double FVolumen0;
	double FVolumen;
	nmVolumenDeposito FVolumenDeposito;
	nmTipoDeposito FTipoDeposito;
	int *FUnion;
	bool asgFUnion;
	bool asgFNumUniones;
	int *FNumConductos;
	bool asgFNumConductos;
	int *FNumNodos;
	bool asgFNumNodos;
	int *FSentidoFlujo;
	bool asgFSentidoFlujo;

	int FNUnionesED;
	int FNumeroUnionesED;

	// La crea pedro. la otra desaparecera
	std::vector<int> FUnionED;
	std::vector<int> FSentidoUED;

	bool FHayCompresor;
	int FNumeroCompresor;
	int FCompresorSentido;

	// Conservative coupling: what each pipe-to-plenum union (and a discharging compressor) moved this update,
	// taken from the pipe's own boundary-face flux (TTubo::TakeEndFlux), signed positive INTO this element.
	std::vector<bool> FUnionFlujoCara;		//!< union i uses the pipe face flux (else legacy BC flow)
	std::vector<double> FUnionMasa;			//!< mass into the element through union i (kg)
	std::vector<double> FUnionH0;			//!< stagnation enthalpy of that mass (J/kg)
	std::vector<std::vector<double> > FUnionMasaEspecie;	//!< species masses into the element (kg)
	bool FCompFlujoCara;					//!< the compressor term uses its outlet-pipe face flux
	double FCompMasa;						//!< mass into the element through the compressor (kg)
	double FCompH0;							//!< stagnation enthalpy of that mass (J/kg)
	std::vector<double> FCompMasaEspecie;	//!< species masses into the element through the compressor (kg)

	// Variables termodinamicas

	double FTemperature;
	double FPressure;
	double FPresionIsen;
	double FDensidad;
	double FAsonido;
	double FMasa0;
	double FMasa;

	double FHeatPower;

	// Variables temporales

	double FTime;
	// bool FUpdated;

	// Salida grafica

	int FNumResMed;

	stResInstantDep FResInstantDep;
	stResMediosDep FResMediosDep;

	// Sensores

	stSensoresDep FSensorDep;

	// Compressor acoplado (si hay)

	TCompresor *FCompresor;

	nmSentidoFlujo FSentidoGasto;

	// Estudio de Estabilidad

	// bool FEstudioEstabilidadRealizado;
	// bool FCalculadoPaso;

	// Variables del transporte de especies.
	int FNumeroEspecies;
	nmTipoCalculoEspecies FCalculoEspecies;
	nmCalculoGamma FCalculoGamma;
	double *FFraccionMasicaEspecie;
	double *FMasaEspecie;
	double FCpMezcla;
	double FCvMezcla;
	double FRMezcla;
	double FGamma;
	double FGamma1;
	double FGamma2;
	double FGamma3;
	double FGamma5;
	double FGamma6;
	bool FHayEGR;
	int FIntEGR;
	double FFraccionMasicaEspecieFuel; // Para calculo completo de especies

	// ---------------------------------------------------------------------------
	// FUNCIONES PRIVADAS
	// ---------------------------------------------------------------------------

	// void PutNUniones(int value);

	double EntalpiaEntrada(double ASonidoE, double VelocidadE, double MasaE, double ASonidoD, double MasaD, double Gamma);

	/*! Same energy term as EntalpiaEntrada, from the stagnation enthalpy h0 (J/kg) of the exchanged mass MasaE
	 (signed, + into the element); ASonidoD is the element's dimensionless speed of sound. Applied to both flow
	 directions it conserves energy exactly with the pipe's face energy flux. */
	double EntalpiaEntradaH0(double h0, double MasaE, double ASonidoD, double MasaD, double Gamma);

	//! Fill the FUnion*/FComp* caches by draining the connected pipes' boundary-face accumulators (once per update).
	void TomaFlujosUniones();

	// void PutCalculadoPaso(bool valor);

  public:
	// ---------------------------------------------------------------------------
	// VARIABLES PUBLICAS
	// ---------------------------------------------------------------------------

	int FNumeroUniones;

	void PutNUniones(int value);

	int getNUniones() {
		return FNumeroUniones;
	}

	// Numero de uniones total del deposito.
	int getNumeroDeposito() {
		return FNumeroDeposito;
	}
	;

	int getNUnionesED() {
		return FNumeroUnionesED;
	}
	;

	double getMasa() {
		return FMasa;
	}
	;

	double getPressure() {
		return FPressure;
	}
	;

	double getPresionIs() {
		return FPresionIsen;
	}
	;

	double getSpeedsound() {
		return FAsonido;
	}
	;

	double getTemperature() {
		return FTemperature;
	}
	;

	TCondicionContorno* GetCCUnionEntreDep(int i);

	TCondicionContorno* GetCCDeposito(int i);
	double GetFraccionMasicaEspecie(int i);

	// ---------------------------------------------------------------------------
	// ---------------------------------------------------------------------------
	double getR() {
		return FRMezcla;
	}
	;

	double getGamma() {
		return FGamma;
	}
	;

	double getTiempo() {
		return FTime;
	}
	;

	void PutUnion(int indice, int valor); // Asigna los numeros de las uniones a depositos del deposito

	void PutNumConductos(int indice, int valor); // Asigna los numeros de los tubos que descargan al deposito

	void PutNumNodos(int indice, int valor); // Asigna los numeros de los nodos que descargan al deposito

	void PutSentidoFlujo(int indice, int valor); // Asigna el sentido del flujo de cada union

	nmTipoDeposito getTipoDeposito() {
		return FTipoDeposito;
	}
	;

	bool FUpdated;

	bool getUpdated() {
		return FUpdated;
	}

	void PutUpdate(bool valor) {
		FUpdated = valor;
	}

	bool FEstudioEstabilidadRealizado;

	bool getEstudioEstabilidadRealizado() {
		return FEstudioEstabilidadRealizado;
	}

	void PutRealizado(bool valor) {
		FEstudioEstabilidadRealizado = valor;
	}

	void PutCalculadoPaso(bool valor) {
		FCalculadoPaso = valor;
	}

	bool FCalculadoPaso;

	bool getCalculadoPaso() {
		return FCalculadoPaso;
	}

	// ---------------------------------------------------------------------------
	// FUNCIONES PUBLICAS
	// ---------------------------------------------------------------------------

	TDeposito(int i, nmTipoDeposito TipoDeposito, nmTipoCalculoEspecies SpeciesModel, int numeroespecies,
			  nmCalculoGamma GammaCalculation, bool ThereIsEGR);

	virtual ~TDeposito() = 0;

	virtual void ActualizaPropiedades(double TimeCalculo) = 0;

	void ReadInstantaneousResultsDep(const char *FileWAM, fpos_t &filepos);

	void HeaderInstantaneousResultsDep(std::stringstream& insoutput, stEspecies *DatosEspecies);

	void ImprimeResultadosInstantaneosDep(std::stringstream& insoutput);

	void ResultadosInstantaneosDep();

	void ReadAverageResultsDep(const char *FileWAM, fpos_t &filepos);

	void HeaderAverageResultsDep(std::stringstream& medoutput, stEspecies *DatosEspecies);

	void ImprimeResultadosMediosDep(std::stringstream& medoutput);

	void AcumulaResultadosMedios(double Actual);

	void ResultadosMediosDep();

	double CriterioEstabilidad(double TMinimo);

	void SalidaGeneralDep(stEspecies *DatosEspecies);

	void CreaUnionED(int numero, int sentido);

	void AsignaCompresor(TCompresor *Compressor, int sentido);

	void AsignacionCC(TCondicionContorno **BC, int numCC);

	void LeeDatosGeneralesDepositos(const char *FileWAM, fpos_t &filepos);

	void ActualizaTiempo(double TiempoActual);

	virtual void UpdateProperties0DModel(double Time) {
	}
	;

	void PutHeatPower(double Power) {
		FHeatPower = Power;
	}
	;

};

// ---------------------------------------------------------------------------
#endif

