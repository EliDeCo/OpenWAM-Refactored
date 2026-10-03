/*--------------------------------------------------------------------------------*\
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


 \*--------------------------------------------------------------------------------*/

//---------------------------------------------------------------------------
#pragma hdrstop

#include "TController.h"
#include "TPIDController.h"

//---------------------------------------------------------------------------

//stGainInput::stGainInput(){}

TController::TController(nmControlMethod meth, int i) {
	// Value-initialise every struct member first: several flags/accumulators were never set.
	FResMediosCtrl = stResMediosCtrl();
	FResInstantCtrl = stResInstantCtrl();
	FControl = meth;
	FControllerID = i + 1;

	// Zero every flag/accumulator: derived controllers only switch on the results they read, and an
	// uninitialised flag randomly added extra result columns and accumulations.
	FResMediosCtrl = stResMediosCtrl();
	FResInstantCtrl = stResInstantCtrl();
	FResMediosCtrl.Output = false;
	FResMediosCtrl.Error = false;

	FResInstantCtrl.Output = false;
	FResInstantCtrl.Error = false;

}

TController::~TController() {

}

#pragma package(smart_init)
