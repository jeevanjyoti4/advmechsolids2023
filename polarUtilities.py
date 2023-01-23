import sympy as sym
r = sym.symbols('r')
theta = sym.symbols('theta')

Q = sym.Matrix([[sym.cos(theta), sym.sin(theta)],[-sym.sin(theta), sym.cos(theta)]])


delr_delx = sym.cos(theta)
deltheta_delx = -sym.sin(theta)/r
delr_dely = sym.sin(theta)
deltheta_dely = sym.cos(theta)/r


def del_delx(f):
    return (sym.diff(f,r)*delr_delx + sym.diff(f,theta)*deltheta_delx).simplify()

def del_dely(f):
    return (sym.diff(f,r)*delr_dely + sym.diff(f,theta)*deltheta_dely).simplify()


def del2_delx2(f):
    return del_delx(del_delx(f)).simplify()

def del2_dely2(f):
    return del_dely(del_dely(f)).simplify()

def del2_delxdely(f):
    return del_delx(del_dely(f)).simplify()


def polarLaplacian(f):
    return (del2_delx2(f) + del2_dely2(f)).simplify()


def polarBiharmonic(f):
    return polarLaplacian(polarLaplacian(f)).simplify()


def sigma_xx(f):
    return del2_dely2(f).simplify().expand()

def sigma_yy(f):
    return del2_delx2(f).simplify().expand()

def sigma_xy(f):
    return -del2_delxdely(f).simplify().expand()


def sigma_Cart(f):
    return sym.Matrix([[sigma_xx(f), sigma_xy(f)],[sigma_xy(f), sigma_yy(f)]])


def sigma_Polar(f):
    return Q*sigma_Cart(f)*(Q.T)


def sigma_rr(f):
    return (sigma_Polar(f)[0,0]).simplify().expand()

def sigma_tt(f):
    return (sigma_Polar(f)[1,1]).simplify().expand()

def sigma_rt(f):
    return (sigma_Polar(f)[0,1]).simplify().expand()

