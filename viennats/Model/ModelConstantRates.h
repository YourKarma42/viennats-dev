#ifndef MODELCONSTANTRATES_H_
#define MODELCONSTANTRATES_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <stack>
#include <vector>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include <iostream>

using namespace std;

namespace model {

	class ConstantRates {

        double StartDirection[3];
        std::vector<double> isotropic_rates;
        std::vector<double> directional_rates;
        std::vector<double> constant_rates;
        std::vector<double> vector_rates;


	public:

        bool MaskLayer;
        bool DirEmpty;
        static const bool SpatiallyEqualDistributedFlux=true;
	    static const bool ReemissionIsMaterialDependent=false;
	    bool CalculateConnectivities;
	    bool CalculateVisibilities;
	    bool CalculateNormalVectors;
	    bool IncludeVectorRates;

        static const bool DropletDistribution=false;
        static const unsigned int NumberOfDroplets=0;

        class DropletType {
        public:
            double Velocity[3];
            double Position[3];
            double Radius;
            double Charge;
        };

		class ParticleType {
		public:
			double Direction[3];
			double Flux;
		};

//		ConstantRates(const std::string & Parameters) : CalculateConnectivities(true), CalculateVisibilities(true),CalculateNormalVectors(true) {
//		    using namespace boost::spirit::classic;
//
//            bool b = parse(
//                    Parameters.begin(),
//                    Parameters.end(),
//                    *(
//                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
//                            (str_p("constant_rates")  >> '='  >>  '{' >> (real_p[push_back_a(constant_rates)] % ',') >> '}'  >> ';') |
//                            (str_p("isotropic_rates")  >> '='  >>  '{' >> (real_p[push_back_a(isotropic_rates)] % ',') >> '}'  >> ';') |
//                            (str_p("directional_rates")  >> '='  >>  '{' >> (real_p[push_back_a(directional_rates)] % ',') >> '}'  >> ';')
//                    ),
//                    space_p | comment_p("//") | comment_p("/*", "*/")).full;
//
//            if (isotropic_rates.size()==0) CalculateConnectivities=false;
//            if (directional_rates.size()==0) {
//                CalculateVisibilities=false;
//                CalculateNormalVectors=false;
//            }
//
//            if (!b) msg::print_error("Failed interpreting process parameters!");
//
//
//		}

		ConstantRates(const std::string & Parameters, bool masked=false) : CalculateConnectivities(true), CalculateVisibilities(true),CalculateNormalVectors(true),IncludeVectorRates(true) {
		    using namespace boost::spirit::classic;

		    bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("constant_rates")  >> '='  >>  '{' >> (real_p[push_back_a(constant_rates)] % ',') >> '}'  >> ';') |
                            (str_p("isotropic_rates")  >> '='  >>  '{' >> (real_p[push_back_a(isotropic_rates)] % ',') >> '}'  >> ';') |
                            (str_p("directional_rates")  >> '='  >>  '{' >> (real_p[push_back_a(directional_rates)] % ',') >> '}'  >> ';') |
                            (str_p("vector_rates")  >> '='  >>  '{' >> (real_p[push_back_a(vector_rates)] % ',') >> '}'  >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (isotropic_rates.size()==0) CalculateConnectivities=false;

            if (directional_rates.size()==0) {
                CalculateVisibilities=false;
                if (vector_rates.size()==0) {
                	CalculateNormalVectors=false;
                }
            }
            //if (CalculateNormalVectors) std::cout << "dr.size: " << directional_rates.size() << ", vr.size: " << vector_rates.size() << std::endl;

            if (vector_rates.size()==0) IncludeVectorRates=false;

            DirEmpty=true;
            for(int i=0;i<3;i++) if(StartDirection[i]!=0) DirEmpty=false;
            if (!DirEmpty) {
            	double start_norm=0;
            	for (int i=0;i<3;++i) start_norm+=StartDirection[i]*StartDirection[i];
            	start_norm=std::sqrt(start_norm);
            	for (int i=0;i<3;++i) StartDirection[i]/=start_norm;
            }// else DirEmpty=0;

            MaskLayer=masked;
            if (!b) msg::print_error("Failed interpreting process parameters!");

            // ---------------------------
            //CalculateNormalVectors=true;
            // ---------------------------

		}

		static const int CoverageStorageSize=0;
		static const int RatesStorageSize=0;
		static const unsigned int NumberOfParticleTypes=0;
		static const unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

		template<class VecType>
		void CalculateVelocity(
				double &Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible) const {

		    double isotropic_rate=(Material < static_cast<int>(isotropic_rates.size()))?isotropic_rates[Material]:0;
		    double directional_rate=(Material < static_cast<int>(directional_rates.size()))?directional_rates[Material]:0;
		    double constant_rate=(Material < static_cast<int>(constant_rates.size()))?constant_rates[Material]:0;

		    Velocity=constant_rate;

//		    if (Velocity >= 0) std::cout << "Velocity 1: " << Velocity << " ---> 1";
	    	//std::cout << "Material: " << Material << std::endl;
		    //std::cout << "constant_rates: " << constant_rates[Material] << std::endl;

		    if (connected) Velocity+=isotropic_rate;

		    if ((visible) && (directional_rate!=0)) {
		        double dot=0.;
		        for (int i=0;i<3;++i) dot-=StartDirection[i]*NormalVector[i];
		        Velocity+=directional_rate*std::max(0.,dot);
		    }
//		    if (Velocity >= 0) std::cout << "Velocity 2: " << Velocity << "\n";

//		    if (IncludeVectorRates) {
//			    double dot2=0.;
//   		        for (int i=0;i<3;++i) {
//   		        	dot2+=vector_rates[i]*NormalVector[i];
//   		        }
//   		        	Velocity+=(Material==0)?dot2:0;
//		    }

            // ---------------------------
//	    	if (abs(NormalVector[1])<0.005){
//	    		Velocity*=10;
//		    }
            // ---------------------------

		}

		template<class VecType>
		void CalculateVectorVelocity(
				double *Velocity,
				double location,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible,
				double oldVel,
				bool Mask) const {

//		    double vector_rate=(Material < static_cast<int>(vector_rates.size()))?vector_rates[Material]:0;
//
//
//			if (IncludeVectorRates){
//				if (MaskLayer) {
//				    double loc=0.05*location;
//			    	loc=(loc>7.5)?15-loc:loc;
//				    loc*=0.01*loc;
//
//for (int i=0;i<3;++i) Velocity[i]=(DirEmpty)?vector_rate*loc:StartDirection[i]*vector_rate*loc;
//				} else {
//for (int i=0;i<3;++i) Velocity[i]=(DirEmpty)?vector_rate:StartDirection[i]*vector_rate;
//				}
//			}

//			double directional_rate=(Material < static_cast<int>(directional_rates.size()))?directional_rates[Material]:0;
//
//			if ((visible) && (directional_rate!=0)) {
//				double dot=0.;
//				for (int i=0;i<3;++i) dot-=StartDirection[i]*NormalVector[i];
//				if (dot>=0.) {
//					for (int i=0;i<3;++i) Velocity[i]=-StartDirection[i]*directional_rate;
//				}
//			}
		}

		static void UpdateCoverage(double *Coverages, const double *Rates) {}

        template <class DropletType, class ParameterType, class PartitionType>
        void DropletGeneration(DropletType& d, double* Position, const ParameterType& Parameter, const PartitionType& Partition) const {}

        template <class PT, class DT, class ParameterType, class PartitionType>
        void ParticleGeneration(PT& p, DT& d, int ParticleType, double ProcessTime, double* Position, const ParameterType& Parameter, const PartitionType& Partition) const {}

//        template <class PT> static void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) {}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {}


		template <class PT, class VecType> static void ParticleReflexion(
							const PT& p,
							std::stack<PT>& particle_stack,
							const VecType& NormalVector,
							const double* Coverages,
							int Material) {}
	};

	const unsigned int ConstantRates::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}

#endif /*MODELCONSTANTRATES_H_*/
