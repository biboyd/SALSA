from os.path import join, exists

import yt
import salsa
import numpy as np

from pytest import fixture

@fixture(scope="module")
def tempdir(tmp_path_factory):
    return tmp_path_factory.mktemp("test_plot_enzo_rays")

@fixture(scope="module")
def extractor(tempdir):
    np.random.seed(567)

    ds = yt.load_sample("IsolatedGalaxy")
    salsa.generate_lrays(ds, 
                         ray_directory=tempdir,
                         center=[0.5, 0.5, 0.5],
                         ion_list=["H I"],
                         n_rays=1,
                         min_impact_param=0,
                         max_impact_param=50)
    
    absext = salsa.SPICEAbsorberExtractor(ds)
    absext.load_ray(join(tempdir, "ray0.h5"))
    absext.get_current_absorbers()
    
    yield absext

    ds.close()

def test_AbsorberPlotter_init(extractor):
    pltr = salsa.AbsorberPlotter(extractor)

def test_AbsorberPlotter_create_slice(extractor, tempdir):
    pltr = salsa.AbsorberPlotter(extractor)
    pltr.create_slice()
    img_path = tempdir / "test_slice.png"
    pltr.slice.save(img_path)
    assert exists(img_path)

def test_AbsorberPlotter_plot_vel_space(extractor):
    pltr = salsa.AbsorberPlotter(extractor)
    vel, flx = pltr.plot_vel_space()

    vel2 = np.arange(-1500, 1510, 10)
    assert np.allclose(vel.value, vel2)

    flx_part = np.array([
       9.99235756e-001, 9.96162370e-001, 9.93731350e-001, 9.71636817e-001,
       4.08938925e-001, 8.97053914e-008, 4.32747229e-052, 5.23372281e-155,
       2.05445475e-212, 3.47363683e-180, 2.61460912e-149, 4.71140332e-143,
       5.50230312e-175, 1.18264905e-202, 2.50923135e-163, 1.65565809e-094,
       8.08243928e-073, 1.41458067e-131, 1.34415203e-216, 1.25754817e-243,
       6.10014581e-226])
    assert np.unique(flx[:135]) == 1
    assert np.allclose(flx[135:156], flx_part)
    assert np.unique(flx[156:242] == 0)
    assert np.unique(flx[242:] == 1)

def test_AbsorberPlotter_plot_lambda_space(extractor):
    pltr = salsa.AbsorberPlotter(extractor)
    wave, flx = pltr.plot_lambda_space()

    wave2 = np.arange(1200.7, 1230.8, .1)
    np.allclose(wave.value, wave2)

    flx_part = np.array([
       9.99005211e-001, 9.94979861e-001, 3.71879654e-003, 3.72195551e-117,
       2.82429723e-178, 2.41024530e-163, 1.50663817e-154, 2.22166394e-108,
       3.99422265e-234, 1.52292551e-134, 5.41625901e-007, 9.86646095e-001,
       9.96252859e-001, 9.98978098e-001])
    assert np.unique(flx[:143]) == 1
    assert np.allclose(flx[143:157], flx_part)
    assert np.unique(flx[157:]) == 1

def test_AbsorberPlotter_plot_multiplot(extractor, tempdir):
    """
    Combines the plotting functions tested above,
    supplying an `ax` argument to plot_vel_space and
    plot_lambda_space. Also tests plot_num_density,
    which doesn't return anything but requires
    1 or 2 axis arguments.
    """
    pltr = salsa.AbsorberPlotter(extractor)
    img_path = tempdir/"test_multiplot.png"
    pltr.plot_multiplot(outfname = img_path)
    assert exists(img_path)
