from numpy import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter as color_conv

num_bins = 40
prefix   = 'smits_{0}'.format(num_bins)

# Plotting colors.
solarized = {
    'base03':  '#002b36',
    'base02':  '#073642',
    'base01':  '#586e75',
    'base00':  '#657b83',
    'base0':   '#839496',
    'base1':   '#93a1a1',
    'base2':   '#eee8d5',
    'base3':   '#fdf6e3',
    'white':   '#fdf6e3',
    'yellow':  '#b58900',
    'orange':  '#cb4b16',
    'red':     '#dc322f',
    'magenta': '#d33682',
    'violet':  '#6c71c4',
    'blue':    '#268bd2',
    'cyan':    '#2aa198',
    'green':   '#859900',
}

# Read the cmf matrix from the given CIE file.
# The file format must be wavelength,x,y,z
def read_cmf(filename, minWave, maxWave):
    import re

    expr = re.compile('^(\d+),([\dE.+-]+),([\dE.+-]+),([\dE.+-]+)\s*$')
    f = open(filename, "r")
    wave = []
    x    = []
    y    = []
    z    = []
    nx = 0
    ny = 0
    nz = 0
    for l in f:
        match = expr.search(l)
        if match:
            g = match.groups()
            w_ = int(g[0])
            if w_ < minWave or w_ > maxWave:
                continue
            x_ = float(g[1])
            y_ = float(g[2])
            z_ = float(g[3])
            wave.append(w_)
            x.append(x_)
            nx += x_
            y.append(y_)
            ny += y_
            z.append(z_)
            nz += z_

    deltaL = wave[1]-wave[0]

    # Normalize functions!
    n = max([nx, ny, nz]) * deltaL
    for i in range(len(wave)):
        x[i] = x[i] / n
        y[i] = y[i] / n
        z[i] = z[i] / n

    return (wave, x, y, z)

# Resample CMF to fewer bins.
def resample(wave, x, y, z, n):
    w_step = wave[1]-wave[0]
    w_min = wave[0]
    w_max = w_min + len(wave) * w_step

    wave_ = linspace(w_min, w_max, n)
    x_ = interp(wave_, wave, x)
    y_ = interp(wave_, wave, y)
    z_ = interp(wave_, wave, z)

    return (wave_, x_, y_, z_)

cmf_orig = [[],[],[]]
cmf      = [[],[],[]]
wave_orig, cmf_orig[0], cmf_orig[1], cmf_orig[2] = read_cmf('ciexyz31.csv', 380, 780)
wave, cmf[0], cmf[1], cmf[2] = resample(wave_orig, cmf_orig[0], cmf_orig[1], cmf_orig[2], num_bins)

# Plot all CMF.

functions = [
    { 'name': 'x', 'color': 'red' },
    { 'name': 'y', 'color': 'green' },
    { 'name': 'z', 'color': 'blue' },
]

plt.figure(1, figsize=(24, 4.5))

for i in range(len(functions)):
    plt.subplot(1, len(functions), i+1)
    f = functions[i]
    plt.fill_between(wave_orig, cmf_orig[i], color=color_conv.to_rgba(solarized[f['color']], 0.5))

    line, = plt.plot(wave_orig, cmf_orig[i])
    plt.setp(line, color=solarized[f['color']], linewidth=2)

    line, = plt.plot(wave, cmf[i])
    plt.setp(line, color=solarized['base3'], linewidth=2)

    plt.gca().set_axis_bgcolor(solarized['base02'])
    plt.gca().grid(True, color=solarized['base01'])

    plt.xlabel('$\lambda$')
    plt.ylabel('$\\bar{' + f['name'] + '}(\lambda)$')
    plt.ylim([min(map(min, cmf)), max(map(max, cmf))*1.1])

plt.savefig(prefix + '_cmf.svg')
plt.savefig(prefix + '_cmf.png')

# Generate a difference vector (abs(p0-p1), abs(p1-p2), ...)
# adjacent values in a vector.
def roughness(p):
    r = 0
    for i in range(len(p)-1):
        r += (p[i]-p[i+1])**2
    return sqrt(r)

# Size of a wavelength bin.
wave_step = wave[1]-wave[0]

# Matrix to convert from spectrum to xyz.
xyz_from_spectrum = wave_step * matrix(cmf)

# Matrix to convert from xyz to srgb.
srgb_from_xyz = matrix([
    [3.2406, -1.5372, -0.4986],
    [-0.9689,  1.8758,  0.0415],
    [0.0557,  -0.204,   1.0570]
])

# Colors to be converted.
colors = {
    'white':   matrix('1; 1; 1'),
    'cyan':    matrix('0; 1; 1'),
    'magenta': matrix('1; 0; 1'),
    'yellow':  matrix('1; 1; 0'),
    'red':     matrix('1; 0; 0'),
    'green':   matrix('0; 1; 0'),
    'blue':    matrix('0; 0; 1'),
}

spectra = {}

# Compute optimal fit for all the colors above.
for name, color in colors.items():
    x0 = matrix([0] * len(wave)).transpose()

    cstr = (
        {
            'type': 'eq',
            'fun': lambda s: (srgb_from_xyz * xyz_from_spectrum * s.reshape(num_bins, 1) - color).transpose().tolist()[0]
        },
    )
    bnds = tuple([(0, 1) for i in range(len(wave))])

    res = minimize(roughness, x0, method='SLSQP', constraints=cstr, bounds=bnds)

    spectra[name] = res.x

plt.figure(2, figsize=(8, 6))

def print_cpp_array(name, l):
    out = 'std::array<float, {0}> {1} = \n\t {{{{'.format(len(l), name)
    out += ', '.join(map(str, l))
    out += '}};'
    print(out)

# Output spectra.
print_cpp_array('wavelengths', wave)
print_cpp_array('cmf_x', cmf[0])
print_cpp_array('cmf_y', cmf[1])
print_cpp_array('cmf_z', cmf[2])
for name, color in colors.items():
    print_cpp_array(name, spectra[name])

# Display spectra.
for name, color in colors.items():
    line, = plt.plot(wave, spectra[name])
    plt.setp(line, color=solarized[name], linewidth=2)

    plt.xlabel('$\lambda$')
    plt.ylabel('Power')
    plt.gca().set_axis_bgcolor(solarized['base02'])
    plt.gca().grid(True, color=solarized['base01'])

plt.savefig(prefix + '_spectra.png')
plt.savefig(prefix + '_spectra.svg')
