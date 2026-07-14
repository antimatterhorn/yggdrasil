import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from yggdrasil import *
from EOS import *
from MieGruneisenMaterials import *
from TillotsonMaterials import *
from typing import Any, Callable

EOSFunc = Callable[[Any], Any]

def IdealGasFunc(constants):
    return IdealGasEOS(1.4, constants)

def HelmholtzFunc(constants):
    return HelmholtzEOS("genHelm.dat", constants)

def MieGruneisenFunc(constants):
    params = MieGruneisenMaterial("iron")  # still in CGS
    return MieGruneisenEOS(constants=constants, **params)

def PolytropicFunc(constants):
    return PolytropicEOS(1.0, 1.4, constants)

def TillotsonFunc(constants):
    params = TillotsonMaterial("iron")  # still in CGS
    return TillotsonEOS(constants=constants, **params)

def IsothermalFunc(constants):
    return IsothermalEOS(5.0, constants)

getEOS: dict[str, EOSFunc] = {
    "IdealGasEOS": IdealGasFunc,
    "HelmholtzEOS": HelmholtzFunc,
    "MieGruneisenEOS": MieGruneisenFunc,
    "PolytropicEOS": PolytropicFunc,
    "TillotsonEOS": TillotsonFunc,
    "IsothermalEOS": IsothermalFunc,
}

if __name__ == "__main__":
    commandLine = CommandLineArguments(nrho = 50,
                                       nu = 50,
                                       minrho = 1e-2,
                                       maxrho = 1e6,
                                       minT   = 200,
                                       maxT   = 1e8,
                                       eos = "IdealGasEOS")

    assert eos in getEOS, "unknown eos '%s', must be one of %s" % (eos, sorted(getEOS.keys()))

    constants = CGS()
    minu = constants.kB * minT / (0.6 * constants.protonMass)
    maxu = constants.kB * maxT / (0.6 * constants.protonMass)
    # Log-space grid
    log_rho = np.linspace(np.log10(minrho), np.log10(maxrho), nrho)
    log_u   = np.linspace(np.log10(minu),   np.log10(maxu),   nu)
    RHO, U  = np.meshgrid(10**log_rho, 10**log_u, indexing='ij')

    # Flattened input for NodeList and Fields
    N = nrho * nu

    eos = getEOS[eos](constants)

    print(eos)
    print(type(eos))

    print(eos.name())

    rho_field = FieldofDouble("rho")
    u_field   = FieldofDouble("internalEnergy")
    p_field   = FieldofDouble("pressure")

    for idx in range(N):
        i = idx // nu
        j = idx % nu
        rho_val = RHO[i, j]
        u_val   = U[i, j]
        rho_field.addValue(rho_val)
        u_field.addValue(u_val)
        p_field.addValue(0) # fill in with eos call


    # EOS computation
    eos.setPressure(p_field, rho_field, u_field)

    # Reformat output to 2D array for plotting
    P = np.array([p_field[i] for i in range(N)]).reshape((nrho, nu))

    # 3D surface plot
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(np.log10(RHO), np.log10(U), np.log10(P),
                    cmap=cm.viridis, edgecolor='none')
    ax.set_xlabel("log10(Rho)")
    ax.set_ylabel("log10(U)")
    ax.set_zlabel("log10(P)")
    ax.set_title(eos.name())

    plt.tight_layout()
    plt.show()

