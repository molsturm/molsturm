include(ExternalData)
set(ExternalData_URL_TEMPLATES "https://get.molsturm.org/data_files/molsturm/%(algo)/%(hash)")

ExternalData_Expand_Arguments(
molsturm-download-testdata
DOWNLOADED

DATA{./interface_python/testdata/be_ccpvtz.hf.hdf5}
DATA{./interface_python/testdata/be_ccpvqz.hf.hdf5}
DATA{./interface_python/testdata/be_ccpvdz.hf.hdf5}
)
ExternalData_Add_Target(molsturm-download-testdata)
