# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# compile CXX with /Library/Developer/CommandLineTools/usr/bin/c++
CXX_FLAGS = -W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -Wno-variadic-macros -Wshadow -pipe -Qunused-arguments -DGL_SILENCE_DEPRECATION -stdlib=libc++ -DG4USE_STD11 -std=c++11 -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fPIE -fstack-protector-strong -O2 -pipe -stdlib=libc++ -fvisibility-inlines-hidden  -fmessage-length=0 -I/opt/anaconda3/envs/root6/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/root_1607274702958/work=/usr/local/src/conda/root_base-6.22.6 -fdebug-prefix-map=/opt/anaconda3/envs/root6=/usr/local/src/conda-prefix -std=c++17 -m64 -pipe -fsigned-char -fno-common -Qunused-arguments -pthread -stdlib=libc++ -O3 -DNDEBUG -fno-trapping-math -ftree-vectorize -fno-math-errno -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk   -fPIC -std=c++11

CXX_DEFINES = -DG4INTY_USE_QT -DG4INTY_USE_XT -DG4UI_USE -DG4UI_USE_QT -DG4UI_USE_TCSH -DG4VERBOSE -DG4VIS_USE -DG4VIS_USE_OPENGL -DG4VIS_USE_OPENGLQT -DG4VIS_USE_OPENGLX -DG4_STORE_TRAJECTORY -DQT_CORE_LIB -DQT_GUI_LIB -DQT_NO_DEBUG -DQT_OPENGL_LIB -DQT_PRINTSUPPORT_LIB -DQT_WIDGETS_LIB

CXX_INCLUDES = -I/opt/anaconda3/envs/root6/include -I/Users/vincentgousy-leblanc/Desktop/Research/Geant4/Geant4_PTF_2/Geant4_PTF/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/g4tools/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/accumulables/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/csv/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/hntools/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/root/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/analysis/xml/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/digits_hits/detector/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/digits_hits/digits/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/digits_hits/hits/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/digits_hits/scorer/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/digits_hits/utils/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/error_propagation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/event/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/externals/clhep/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/externals/zlib/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/biasing/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/divisions/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/magneticfield/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/navigation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/solids/Boolean/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/solids/CSG/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/solids/specific/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/geometry/volumes/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/global/HEPGeometry/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/global/HEPNumerics/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/global/HEPRandom/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/global/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/graphics_reps/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/intercoms/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/interfaces/GAG/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/interfaces/basic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/interfaces/common/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/materials/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/parameterisations/gflash/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/adjoint/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/bosons/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/hadrons/barions/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/hadrons/ions/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/hadrons/mesons/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/leptons/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/shortlived/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/particles/utils/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/persistency/ascii/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/persistency/mctruth/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/builders/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/decay/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/electromagnetic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/factory/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/gamma_lepto_nuclear/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/hadron_elastic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/hadron_inelastic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/ions/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/limiters/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/constructors/stopping/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/lists/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/physics_lists/util/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/biasing/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/biasing/generic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/biasing/importance/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/cuts/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/decay/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/adjoint/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/dna/processes/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/dna/models/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/dna/utils/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/dna/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/dna/molecules/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/dna/molecules/types/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/highenergy/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/lowenergy/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/muons/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/pii/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/polarisation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/standard/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/utils/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/electromagnetic/xrays/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/cross_sections/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/abla/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/abrasion/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/binary_cascade/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/cascade/cascade/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/coherent_elastic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/ablation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/evaporation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/fermi_breakup/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/fission/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/gem_evaporation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/handler/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/multifragmentation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/photon_evaporation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/de_excitation/util/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/em_dissociation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/fission/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/im_r_matrix/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/inclxx/utils/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/inclxx/incl_physics/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/inclxx/interface/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/gamma_nuclear/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/lend/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/lepto_nuclear/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/particle_hp/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/parton_string/diffraction/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/parton_string/hadronization/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/parton_string/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/parton_string/qgsm/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/pre_equilibrium/exciton_model/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/qmd/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/quasi_elastic/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/radioactive_decay/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/rpg/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/theo_high_energy/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/models/util/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/processes/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/stopping/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/hadronic/util/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/optical/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/solidstate/phonon/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/solidstate/channeling/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/parameterisation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/scoring/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/processes/transportation/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/readout/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/run/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/track/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/tracking/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/FukuiRenderer/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/HepRep/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/RayTracer/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/Tree/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/VRML/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/XXX/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/externals/gl2ps/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/gMocren/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/management/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/modeling/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/source/visualization/OpenGL/include -isystem /Users/vincentgousy-leblanc/Downloads/geant4.10.04.p03/build/source/externals/zlib -isystem /opt/anaconda3/include/qt -isystem /opt/anaconda3/include/qt/QtWidgets -isystem /opt/anaconda3/include/qt/QtGui -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk/System/Library/Frameworks/OpenGL.framework/Headers -isystem /opt/anaconda3/include/qt/QtCore -isystem /opt/anaconda3/./mkspecs/macx-clang -isystem /opt/anaconda3/include/qt/QtPrintSupport -isystem /opt/anaconda3/include/qt/QtOpenGL 

