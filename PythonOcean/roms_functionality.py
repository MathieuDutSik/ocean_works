import math
from math import cos, sin, acos, floor
from math import sqrt, asin, tan, atan, tanh, sinh, pow, cosh
import numpy as np
from collections import namedtuple
import datetime


ARVDtyp = namedtuple(
    "ARVDtyp", "N Vstretching Vtransform Tcline hc theta_s theta_b Cs_r Cs_w sc_r sc_w"
)
VerticalStrat = namedtuple("VerticalStrat", "Hz z_r z_w")


def RemoveFileIfExist(eFile):
    import os

    if os.path.isfile(eFile):
        os.remove(eFile)


def nearx(LON_arr, eVal):
    len = LON_arr.shape[0]
    idx = len - 1
    eDist = abs(LON_arr[idx] - eVal)
    for i in range(len - 1):
        fDist = abs(LON_arr[i] - eVal)
        if fDist < eDist:
            eDist = fDist
            idx = i
    return idx


def uvp_mask(mask_rho):
    [eta_rho, xi_rho] = mask_rho.shape
    mask_v = np.zeros(shape=(eta_rho - 1, xi_rho))
    mask_u = np.zeros(shape=(eta_rho, xi_rho - 1))
    mask_psi = np.zeros(shape=(eta_rho - 1, xi_rho - 1))
    for i in range(eta_rho - 1):
        for j in range(xi_rho):
            mask_v[i, j] = mask_rho[i, j] * mask_rho[i + 1, j]
    for i in range(eta_rho):
        for j in range(xi_rho - 1):
            mask_u[i, j] = mask_rho[i, j] * mask_rho[i, j + 1]
    for i in range(eta_rho - 1):
        for j in range(xi_rho - 1):
            mask_u[i, j] = (
                mask_rho[i, j]
                * mask_rho[i, j + 1]
                * mask_rho[i + 1, j]
                * mask_rho[i + 1, j + 1]
            )
    return [mask_u, mask_v, mask_psi]


def UTIL_EulerVectorBasis(londeg, latdeg, angle):
    lonRad = londeg * (math.pi / 180)
    latRad = latdeg * (math.pi / 180)
    #
    Ur = [cos(lonRad), sin(lonRad), 0]
    Ulon = [-sin(lonRad), cos(lonRad), 0]
    Uz = [0, 0, 1]
    #
    Urho = [
        cos(latRad) * Ur[0] + sin(latRad) * Uz[0],
        cos(latRad) * Ur[1] + sin(latRad) * Uz[1],
        cos(latRad) * Ur[2] + sin(latRad) * Uz[2],
    ]
    Ulat = [
        -sin(latRad) * Ur[0] + cos(latRad) * Uz[0],
        -sin(latRad) * Ur[1] + cos(latRad) * Uz[1],
        -sin(latRad) * Ur[2] + cos(latRad) * Uz[2],
    ]
    #
    Vect1 = Urho
    Vect2 = [
        cos(anglerad) * Ulon[0] + sin(anglerad) * Ulat[0],
        cos(anglerad) * Ulon[1] + sin(anglerad) * Ulat[1],
        cos(anglerad) * Ulon[2] + sin(anglerad) * Ulat[2],
    ]
    Vect3 = [
        -sin(anglerad) * Ulon[0] + cos(anglerad) * Ulat[0],
        -sin(anglerad) * Ulon[1] + cos(anglerad) * Ulat[1],
        -sin(anglerad) * Ulon[2] + cos(anglerad) * Ulat[2],
    ]
    return [Vect1, Vect2, Vect3]


def UTIL_XYZFromLonLat(lonDeg, latDeg):
    lonRad = lonDeg * (math.pi / 180)
    latRad = latDeg * (math.pi / 180)
    #
    x = cos(lonRad) * cos(latRad)
    y = sin(lonRad) * cos(latRad)
    z = sin(latRad)
    return [x, y, z]


def GetAngle(x, y):
    eComplex = x + y * j
    return np.angle(eComplex, deg=False)


def UTIL_GetLonLat(eX, eY, eZ):
    eNorm = sqrt(eX * eX + eY * eY + eZ * eZ)
    x = eX / eNorm
    y = eY / eNorm
    z = eZ / eNorm
    #
    TheLonRad = GetAngle(x, y)
    #
    rNew = sqrt(x * x + y * y)
    TheLatRad = GetAngle(rNew, z)
    #
    TheLon = TheLonRad * (180 / math.pi)
    TheLat = TheLatRad * (180 / math.pi)
    return [TheLon, TheLat]


def GRID_TransformToLonLat(LonNew, LatNew, V1, V2, V3):
    [eta_rho, xi_rho] = LonNew.shape

    xPosNew = np.zeros(shape=(eta_rho, xi_rho))
    yPosNew = np.zeros(shape=(eta_rho, xi_rho))
    zPosNew = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            [xNew, yNew, zNew] = UTIL_XYZFromLonLat(
                LonNew[iEta, iXi], LatNew[iEta, iXi]
            )
            xPosNew[iEta, iXi] = xNew
            yPosNew[iEta, iXi] = yNew
            zPosNew[iEta, iXi] = zNew
    #
    xPos = np.zeros(shape=(eta_rho, xi_rho))
    yPos = np.zeros(shape=(eta_rho, xi_rho))
    zPos = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            x = (
                xPosNew[iEta, iXi] * V1[0]
                + yPosNew[iEta, iXi] * V2[0]
                + zPosNew[iEta, iXi] * V3[0]
            )
            y = (
                xPosNew[iEta, iXi] * V1[1]
                + yPosNew[iEta, iXi] * V2[1]
                + zPosNew[iEta, iXi] * V3[1]
            )
            z = (
                xPosNew[iEta, iXi] * V1[2]
                + yPosNew[iEta, iXi] * V2[2]
                + zPosNew[iEta, iXi] * V3[2]
            )
            xPos[iEta, iXi] = x
            yPos[iEta, iXi] = y
            zPos[iEta, iXi] = z
    #
    LonMat = np.zeros(shape=(eta_rho, xi_rho))
    LatMat = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            [TheLon, TheLat] = UTIL_GetLonLat(
                xPos[iEta, iXi], yPos[iEta, iXi], zPos[iEta, iXi]
            )
            LonMat[iEta, iXi] = TheLon
            LatMat[iEta, iXi] = TheLat
    return [LonMat, LatMat]


def spheric_dist(lon1, lat1, lon2, lat2):
    # determine the distance on the earth between 2 point
    earthradius = 6356750.52
    deg2rad = math.pi / 180
    rad2deg = 180 / math.pi
    #  Determine proper longitudinal shift.
    delta = lon2 - lon1
    l = abs(delta)
    if l > 180:
        l = 360 - l
    #  Convert Decimal degrees to radians.
    beta1 = lat1 * deg2rad
    beta2 = lat2 * deg2rad
    l = l * deg2rad
    #  Calculate S/Bo subformulas.
    st = sqrt(
        pow(sin(l) * cos(beta2), 2)
        + pow(sin(beta2) * cos(beta1) - sin(beta1) * cos(beta2) * cos(l), 2)
    )
    #       Calculate distance from point 1 to point 2
    dist = asin(st) * earthradius
    return dist


def sign(eX):
    if eX == 0:
        return 0
    if eX > 0:
        return 1
    if eX < 0:
        return -1


def get_angle_corr(LON_u, LAT_u):
    #
    # Compute the grid orientation: angle [radians]
    # between XI-axis and the direction to the EAST
    # at RHO-points.
    #
    # lonu longitude of u points
    # latu latitude of u points
    # argu1: spheroid
    #         'clarke66'  Clarke 1866
    #         'iau73'     IAU 1973
    #         'wgs84'     WGS 1984 (default)
    #         'sphere'    Sphere of radius 6371.0 km
    #
    # copied from dist.m of the Oceans toolbox
    #
    spheroid = "wgs84"

    [eta_u, xi_u] = LON_u.shape
    if spheroid == "sph":
        A = 6371000.0
        B = A
        E = sqrt(A * A - B * B) / A
    if spheroid == "cla":
        A = 6378206.4e0
        B = 6356583.8e0
        E = sqrt(A * A - B * B) / A
    if spheroid == "iau":
        A = 6378160.0e0
        B = 6356774.516e0
        E = sqrt(A * A - B * B) / A
    if spheroid == "wgs84":
        A = 6378137.0
        E = 0.081819191
        B = sqrt(A * A - (A * E) * (A * E))
    EPS = E * E / (1 - E * E)

    LONrad_u = LON_u * (math.pi / 180)
    LATrad_u = LAT_u * (math.pi / 180)

    for iEta in range(eta_u):
        for iXi in range(xi_u):
            if LATrad_u[iEta, iXi] == 0:
                LATrad_u[iEta, iXi] = EPS
    eta_rho = eta_u
    xi_rho = xi_u + 1

    PHI1 = np.zeros(shape=(eta_u, xi_u - 1))
    XLAM1 = np.zeros(shape=(eta_u, xi_u - 1))
    PHI2 = np.zeros(shape=(eta_u, xi_u - 1))
    XLAM2 = np.zeros(shape=(eta_u, xi_u - 1))
    azim3 = np.zeros(shape=(eta_u, xi_u - 1))
    for iEta in range(eta_u):
        for iXi in range(xi_u - 1):
            PHI1[iEta, iXi] = LATrad_u[iEta, iXi]
            XLAM1[iEta, iXi] = LONrad_u[iEta, iXi]
            PHI2[iEta, iXi] = LATrad_u[iEta, iXi + 1]
            XLAM2[iEta, iXi] = LONrad_u[iEta, iXi + 1]
    for iEta in range(eta_u):
        for iXi in range(xi_u - 1):
            if PHI1[iEta, iXi] == PHI2[iEta, iXi]:
                PHI2[iEta, iXi] = PHI2[iEta, iXi] + 1e-14
            if XLAM1[iEta, iXi] == XLAM2[iEta, iXi]:
                XLAM2[iEta, iXi] = XLAM2[iEta, iXi] + 1e-14
            ePHI1 = PHI1[iEta, iXi]
            ePHI2 = PHI2[iEta, iXi]
            eXLAM1 = XLAM1[iEta, iXi]
            eXLAM2 = XLAM2[iEta, iXi]
            xnu1 = A / sqrt(1 - pow(E * sin(ePHI1), 2))
            xnu2 = A / sqrt(1 - pow(E * sin(ePHI2), 2))
            TPSI2 = (1 - E * E) * tan(ePHI2) + E * E * xnu1 * sin(ePHI1) / (
                xnu2 * cos(ePHI2)
            )
            #
            DLAM = eXLAM2 - eXLAM1
            CTA12 = (cos(ePHI1) * TPSI2 - sin(ePHI1) * cos(DLAM)) / sin(DLAM)
            azim = atan(1 / CTA12)
            #
            if abs(DLAM) < math.pi:
                DLAM2 = DLAM
            if DLAM > math.pi:
                DLAM2 = DLAM - 2 * math.pi
            if DLAM < -math.pi:
                DLAM2 = DLAM + 2 * math.pi
            #
            if abs(azim) <= math.pi:
                azim2 = azim
            if azim < -math.pi:
                azim2 = azim + 2 * math.pi
            if azim > math.pi:
                azim2 = azim - 2 * math.pi
            #
            alpha = 0
            if sign(azim2) != sign(DLAM2):
                alpha = 1
            azim3[iEta, iXi] = azim2 + math.pi * sign(-azim2) * alpha
    angle = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_u - 1):
            angle[iEta, iXi + 1] = (math.pi / 2) - azim3[iEta, iXi]
    for iEta in range(eta_rho):
        angle[iEta, 0] = angle[iEta, 1]
        angle[iEta, xi_rho - 1] = angle[iEta, xi_rho - 2]
    return angle


def GRID_DirectWriteGridFile_Generic(
    GridFile,
    lon_rho,
    lat_rho,
    lon_u,
    lat_u,
    lon_v,
    lat_v,
    lon_psi,
    lat_psi,
    f,
    angle,
    pm,
    pn,
    dndx,
    dmde,
    xl,
    el,
    mask_rho,
    h,
):
    import netCDF4
    import numpy as np

    [mask_u, mask_v, mask_psi] = uvp_mask(mask_rho)
    [eta_rho, xi_rho] = mask_rho.shape
    [eta_u, xi_u] = mask_u.shape
    [eta_v, xi_v] = mask_v.shape
    [eta_psi, xi_psi] = mask_psi.shape
    #
    dataset = netCDF4.Dataset(GridFile, "w")
    Nc_eta_rho = dataset.createDimension("eta_rho", eta_rho)
    Nc_eta_u = dataset.createDimension("eta_u", eta_u)
    Nc_eta_v = dataset.createDimension("eta_v", eta_v)
    Nc_eta_psi = dataset.createDimension("eta_psi", eta_psi)
    Nc_xi_rho = dataset.createDimension("xi_rho", xi_rho)
    Nc_xi_u = dataset.createDimension("xi_u", xi_u)
    Nc_xi_v = dataset.createDimension("xi_v", xi_v)
    Nc_xi_psi = dataset.createDimension("xi_psi", xi_psi)
    Nc_one = dataset.createDimension("one", 1)
    #
    NC_dmde = dataset.createVariable("dmde", np.float64, ("eta_rho", "xi_rho"))
    NC_dmde.long_name = "eta derivative of inverse metric factor pm"
    NC_dmde.units = "meter"
    NC_dmde[:, :] = dmde
    #
    NC_dndx = dataset.createVariable("dndx", np.float64, ("eta_rho", "xi_rho"))
    NC_dndx.long_name = "xi derivative of inverse metric factor pn"
    NC_dndx.units = "meter"
    NC_dndx[:, :] = dndx
    #
    NC_el = dataset.createVariable("el", np.float64, ("one"))
    NC_el.long_name = "domain length in the ETA-direction"
    NC_el.units = "degrees"
    NC_el[:] = el
    #
    NC_xl = dataset.createVariable("xl", np.float64, ("one"))
    NC_xl.long_name = "domain length in the XI-direction"
    NC_xl.units = "degrees"
    NC_xl[:] = xl
    #
    NC_f = dataset.createVariable("f", np.float64, ("eta_rho", "xi_rho"))
    NC_f.long_name = "Coriolis parameter at RHO-points"
    NC_f.units = "second-1"
    NC_f[:, :] = f
    #
    NC_h = dataset.createVariable("h", np.float64, ("eta_rho", "xi_rho"))
    NC_h.long_name = "bathymetry at RHO-points"
    NC_h.units = "meter"
    NC_h[:, :] = h
    #
    NC_hraw = dataset.createVariable("hraw", np.float64, ("eta_rho", "xi_rho"))
    NC_hraw.long_name = "raw bathymetry at RHO-points"
    NC_hraw.units = "meter"
    NC_hraw[:, :] = h
    #
    NC_angle = dataset.createVariable("angle", np.float64, ("eta_rho", "xi_rho"))
    NC_angle.long_name = "angle between xi axis and east"
    NC_angle.units = "radiant"
    NC_angle[:, :] = angle
    #
    NC_lat_rho = dataset.createVariable("lat_rho", np.float64, ("eta_rho", "xi_rho"))
    NC_lat_rho.long_name = "latitude of RHO-points"
    NC_lat_rho.units = "degree_north"
    NC_lat_rho[:, :] = lat_rho
    #
    NC_lat_psi = dataset.createVariable("lat_psi", np.float64, ("eta_psi", "xi_psi"))
    NC_lat_psi.long_name = "latitude of PSI-points"
    NC_lat_psi.units = "degree_north"
    NC_lat_psi[:, :] = lat_psi
    #
    NC_lat_u = dataset.createVariable("lat_u", np.float64, ("eta_u", "xi_u"))
    NC_lat_u.long_name = "latitude of U-points"
    NC_lat_u.units = "degree_north"
    NC_lat_u[:, :] = lat_u
    #
    NC_lat_v = dataset.createVariable("lat_v", np.float64, ("eta_v", "xi_v"))
    NC_lat_v.long_name = "latitude of V-points"
    NC_lat_v.units = "degree_north"
    NC_lat_v[:, :] = lat_v
    #
    NC_lon_rho = dataset.createVariable("lon_rho", np.float64, ("eta_rho", "xi_rho"))
    NC_lon_rho.long_name = "longitude of RHO-points"
    NC_lon_rho.units = "degree_east"
    NC_lon_rho[:, :] = lon_rho
    #
    NC_lon_psi = dataset.createVariable("lon_psi", np.float64, ("eta_psi", "xi_psi"))
    NC_lon_psi.long_name = "longitude of PSI-points"
    NC_lon_psi.units = "degree_east"
    NC_lon_psi[:, :] = lon_psi
    #
    NC_lon_u = dataset.createVariable("lon_u", np.float64, ("eta_u", "xi_u"))
    NC_lon_u.long_name = "longitude of U-points"
    NC_lon_u.units = "degree_east"
    NC_lon_u[:, :] = lon_u
    #
    NC_lon_v = dataset.createVariable("lon_v", np.float64, ("eta_v", "xi_v"))
    NC_lon_v.long_name = "longitude of V-points"
    NC_lon_v.units = "degree_east"
    NC_lon_v[:, :] = lon_v
    #
    NC_mask_rho = dataset.createVariable("mask_rho", np.float64, ("eta_rho", "xi_rho"))
    NC_mask_rho.long_name = "longitude of RHO-points"
    NC_mask_rho.option_0 = "land"
    NC_mask_rho.option_1 = "water"
    NC_mask_rho[:, :] = mask_rho
    #
    NC_mask_psi = dataset.createVariable("mask_psi", np.float64, ("eta_psi", "xi_psi"))
    NC_mask_psi.long_name = "longitude of PSI-points"
    NC_mask_psi.option_0 = "land"
    NC_mask_psi.option_1 = "water"
    NC_mask_psi[:, :] = mask_psi
    #
    NC_mask_u = dataset.createVariable("mask_u", np.float64, ("eta_u", "xi_u"))
    NC_mask_u.long_name = "longitude of U-points"
    NC_mask_u.option_0 = "land"
    NC_mask_u.option_1 = "water"
    NC_mask_u[:, :] = mask_u
    #
    NC_mask_v = dataset.createVariable("mask_v", np.float64, ("eta_v", "xi_v"))
    NC_mask_v.long_name = "longitude of V-points"
    NC_mask_v.option_0 = "land"
    NC_mask_v.option_1 = "water"
    NC_mask_v[:, :] = mask_v
    #
    NC_pm = dataset.createVariable("pm", np.float64, ("eta_rho", "xi_rho"))
    NC_pm.long_name = "curvilinear coordinate metric in XI"
    NC_pm.units = "meter-1"
    NC_pm[:, :] = pm
    #
    NC_pn = dataset.createVariable("pn", np.float64, ("eta_rho", "xi_rho"))
    NC_pn.long_name = "curvilinear coordinate metric in ETA"
    NC_pn.units = "meter-1"
    NC_pn[:, :] = pn
    #
    NC_spherical = dataset.createVariable("spherical", "S1", ("one"))
    NC_spherical.long_name = "Grid type logical switch"
    NC_spherical.option_T = "spherical"
    NC_spherical[:] = "T"
    #
    dataset.close()
    print("Loeaving GRID_DirectWriteGridFile_Generic")


def GRID_CreateNakedGrid(
    GridFile, lon, lat, angledeg, XdistMeter, EdistMeter, resolMeter
):
    # GRID_CreateNakedGrid(...
    #    GridFile, lon, lat, angledeg, XdistMeter, EdistMeter, resolMeter)
    #
    # This code is adapted from easygrid.
    # It creates a naked ROMS grid with only the geometrical
    # features:
    # ---no bathymetry
    # ---no mask
    #
    # GridFile is the file of the grid.
    # (lon, lat) are the coordinates of the bottom point.
    # angledeg is the angle in degree for the rotation (AT the point
    #    lon, lat)
    #
    # XdistMeter is the distance in XI direction
    # EdistMeter is the distance in ETA direction
    # resolMeter is the horizontal discretization in XI and ETA directions

    anglerad = angledeg * (math.pi / 180)

    [V1, V2, V3] = UTIL_EulerVectorBasis(lon, lat, anglerad)
    el = EdistMeter
    xl = XdistMeter

    rad2deg = 180 / math.pi

    EarthRadius = 6378137
    AngularXDist = XdistMeter / EarthRadius
    AngularEDist = EdistMeter / EarthRadius
    DeltaEta = resolMeter / EarthRadius
    DeltaXi = DeltaEta
    print("AngularXDist=", AngularXDist)
    print("AngularEDist=", AngularEDist)
    print("DeltaEta=", DeltaEta)
    print("DeltaXi=", DeltaXi)
    #
    xi_rho = 1 + math.floor(AngularXDist / DeltaXi)
    eta_rho = 1 + math.floor(AngularEDist / DeltaEta)
    eta_u = eta_rho
    eta_v = eta_rho - 1
    eta_psi = eta_rho - 1
    xi_u = xi_rho - 1
    xi_v = xi_rho
    xi_psi = xi_rho - 1
    #
    Lp = xi_rho
    Mp = eta_rho
    L = Lp - 1
    M = Mp - 1
    Lm = Lp - 2
    Mm = Mp - 2
    #
    LonN_rho = np.zeros(shape=(eta_rho, xi_rho))
    LatN_rho = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            LonN_rho[iEta, iXi] = rad2deg * iXi * DeltaXi
            LatN_rho[iEta, iXi] = rad2deg * iEta * DeltaEta

    [lon_rho, lat_rho] = GRID_TransformToLonLat(LonN_rho, LatN_rho, V1, V2, V3)
    print("1: eta_rho=", eta_rho, "  xi_rho=", xi_rho)
    print("(1,1)            : lon=", lon_rho[0, 0], " lat=", lat_rho[0, 0])
    print(
        "(1,xi_rho)       : lon=",
        lon_rho[0, xi_rho - 1],
        " lat=",
        lat_rho[0, xi_rho - 1],
    )
    print(
        "(eta_rho,1)      : lon=",
        lon_rho[eta_rho - 1, 0],
        " lat=",
        lat_rho[eta_rho - 1, 0],
    )
    print(
        "(eta_rho,xi_rho) : lon=",
        lon_rho[eta_rho - 1, xi_rho - 1],
        " lat=",
        lat_rho[eta_rho - 1, xi_rho - 1],
    )

    LonN_u = np.zeros(shape=(eta_u, xi_u))
    LatN_u = np.zeros(shape=(eta_u, xi_u))
    for iEta in range(eta_u):
        for iXi in range(xi_u):
            LonN_u[iEta, iXi] = rad2deg * (iXi + 0.5) * DeltaXi
            LatN_u[iEta, iXi] = rad2deg * iEta * DeltaEta
    [lon_u, lat_u] = GRID_TransformToLonLat(LonN_u, LatN_u, V1, V2, V3)

    LonN_v = np.zeros(shape=(eta_v, xi_v))
    LatN_v = np.zeros(shape=(eta_v, xi_v))
    for iEta in range(eta_v):
        for iXi in range(xi_v):
            LonN_v[iEta, iXi] = rad2deg * iXi * DeltaXi
            LatN_v[iEta, iXi] = rad2deg * (iEta + 0.5) * DeltaEta
    [lon_v, lat_v] = GRID_TransformToLonLat(LonN_v, LatN_v, V1, V2, V3)

    LonN_psi = np.zeros(shape=(eta_psi, xi_psi))
    LatN_psi = np.zeros(shape=(eta_psi, xi_psi))
    for iEta in range(eta_psi):
        for iXi in range(xi_psi):
            LonN_psi[iEta, iXi] = rad2deg * (iXi + 0.5) * DeltaXi
            LatN_psi[iEta, iXi] = rad2deg * (iEta + 0.5) * DeltaEta
    [lon_psi, lat_psi] = GRID_TransformToLonLat(LonN_psi, LatN_psi, V1, V2, V3)

    dx = np.zeros(shape=(eta_rho, xi_rho))
    dy = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho - 2):
            dx[iEta, iXi + 1] = spheric_dist(
                lon_u[iEta, iXi],
                lat_u[iEta, iXi],
                lon_u[iEta, iXi + 1],
                lat_u[iEta, iXi + 1],
            )
    for iEta in range(eta_rho):
        dx[iEta, 0] = dx[iEta, 1]
        dx[iEta, xi_rho - 1] = dx[iEta, xi_rho - 2]
    for iEta in range(eta_rho - 2):
        for iXi in range(xi_rho):
            dy[iEta + 1, iXi] = spheric_dist(
                lon_v[iEta, iXi],
                lat_v[iEta, iXi],
                lon_v[iEta + 1, iXi],
                lat_v[iEta + 1, iXi],
            )
    for iXi in range(xi_rho):
        dy[0, iXi] = dy[1, iXi]
        dy[eta_rho - 1, iXi] = dy[eta_rho - 2, iXi]
    pm = np.zeros(shape=(eta_rho, xi_rho))
    pn = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            pm[iEta, iXi] = 1 / dx[iEta, iXi]
            pn[iEta, iXi] = 1 / dy[iEta, iXi]
    #
    dndx = np.zeros(shape=(eta_rho, xi_rho))
    dmde = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho - 2):
        for iXi in range(xi_rho):
            dmde[iEta + 1, iXi] = 0.5 * (1 / pm[iEta + 2, iXi] - 1 / pm[iEta, iXi])
    for iEta in range(eta_rho):
        for iXi in range(xi_rho - 2):
            dndx[iEta, iXi + 1] = 0.5 * (1 / pn[iEta, iXi + 2] - 1 / pn[iEta, iXi])
    #
    h = np.zeros(shape=(eta_rho, xi_rho))
    f = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            f[iEta, iXi] = 2 * (7.29e-5) * sin(lat_rho[iEta, iXi] * (math.pi / 180))
    #
    mask_rho = np.zeros(shape=(eta_rho, xi_rho))
    angleradMat = get_angle_corr(lon_u, lat_u)
    #
    print("lon_rho(min/max)=", lon_rho.min(), " / ", lon_rho.max())
    print("lat_rho(min/max)=", lat_rho.min(), " / ", lat_rho.max())
    #
    GRID_DirectWriteGridFile_Generic(
        GridFile,
        lon_rho,
        lat_rho,
        lon_u,
        lat_u,
        lon_v,
        lat_v,
        lon_psi,
        lat_psi,
        f,
        angleradMat,
        pm,
        pn,
        dndx,
        dmde,
        xl,
        el,
        mask_rho,
        h,
    )


def GeodesicDistance_V2(LonDeg1, LatDeg1, LonDeg2, LatDeg2):
    lon1 = math.pi * LonDeg1 / 180
    lat1 = math.pi * LatDeg1 / 180
    x1 = cos(lon1) * cos(lat1)
    y1 = sin(lon1) * cos(lat1)
    z1 = sin(lat1)
    #
    lon2 = math.pi * LonDeg2 / 180
    lat2 = math.pi * LatDeg2 / 180
    x2 = cos(lon2) * cos(lat2)
    y2 = sin(lon2) * cos(lat2)
    z2 = sin(lat2)
    #
    scalprod = x1 * x2 + y1 * y2 + z1 * z2
    if scalprod > 1:
        return 0
    else:
        return acos(scalprod)


def uvp_lonlat(LON_rho, LAT_rho):
    [eta_rho, xi_rho] = LON_rho.shape
    eta_u = eta_rho
    eta_v = eta_rho - 1
    eta_psi = eta_rho - 1
    xi_u = xi_rho - 1
    xi_v = xi_rho
    xi_psi = xi_rho - 1
    LON_u = np.zeros(shape=(eta_u, xi_u))
    LAT_u = np.zeros(shape=(eta_u, xi_u))
    LON_v = np.zeros(shape=(eta_v, xi_v))
    LAT_v = np.zeros(shape=(eta_v, xi_v))
    LON_psi = np.zeros(shape=(eta_psi, xi_psi))
    LAT_psi = np.zeros(shape=(eta_psi, xi_psi))
    #
    for iEta in range(eta_u):
        for iXi in range(xi_u):
            LON_u[iEta, iXi] = (LON_rho[iEta, iXi] + LON_rho[iEta, iXi + 1]) / 2
            LAT_u[iEta, iXi] = (LAT_rho[iEta, iXi] + LAT_rho[iEta, iXi + 1]) / 2
    #
    for iEta in range(eta_v):
        for iXi in range(xi_v):
            LON_v[iEta, iXi] = (LON_rho[iEta, iXi] + LON_rho[iEta + 1, iXi]) / 2
            LAT_v[iEta, iXi] = (LAT_rho[iEta, iXi] + LAT_rho[iEta + 1, iXi]) / 2
    #
    for iEta in range(eta_psi):
        for iXi in range(xi_psi):
            LON_psi[iEta, iXi] = (LON_u[iEta, iXi] + LON_u[iEta + 1, iXi]) / 2
            LAT_psi[iEta, iXi] = (LAT_u[iEta, iXi] + LAT_u[iEta + 1, iXi]) / 2
    #
    return [LON_u, LON_v, LON_psi, LAT_u, LAT_v, LAT_psi]


def GRID_CreateNakedGridMinMaxLonLat(
    GridFile, MinLon, MaxLon, MinLat, MaxLat, resolMeter
):
    lon = MinLon
    lat = MinLat

    angledeg = 0
    EarthRadius = 6378137

    dist1 = GeodesicDistance_V2(MinLon, MinLat, MinLon, MaxLat)
    dist2 = GeodesicDistance_V2(MinLon, MinLat, MaxLon, MinLat)

    DistMeter_LAT = EarthRadius * dist1
    DistMeter_LON = EarthRadius * dist2
    print("DistMeter_LON=", DistMeter_LON)
    print("DistMeter_LAT=", DistMeter_LAT)

    xi_rho = 1 + round(DistMeter_LON / resolMeter)
    eta_rho = 1 + floor(DistMeter_LAT / resolMeter)
    print("2: eta_rho=", eta_rho, " xi_rho=", xi_rho)
    LON_rho = np.zeros(shape=(eta_rho, xi_rho))
    LAT_rho = np.zeros(shape=(eta_rho, xi_rho))

    el = DistMeter_LAT
    xl = DistMeter_LON
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            LON_rho[iEta, iXi] = MinLon + iXi * (MaxLon - MinLon) / (xi_rho - 1)
            LAT_rho[iEta, iXi] = MinLat + iEta * (MaxLat - MinLat) / (eta_rho - 1)
    print("We have LON_rho, LAT_rho")
    [LON_u, LON_v, LON_psi, LAT_u, LAT_v, LAT_psi] = uvp_lonlat(LON_rho, LAT_rho)
    print("We have LON_uvp, LAT_uvp")
    [eta_u, xi_u] = LON_u.shape
    print("eta_u=", eta_u, " xi_u=", xi_u)
    [eta_v, xi_v] = LON_v.shape
    print("eta_v=", eta_v, " xi_v=", xi_v)
    Lp = xi_rho
    Mp = eta_rho
    L = Lp - 1
    M = Mp - 1
    Lm = Lp - 2
    Mm = Mp - 2

    dx = np.zeros(shape=(eta_rho, xi_rho))
    dy = np.zeros(shape=(eta_rho, xi_rho))
    print("dx, dy, pm, pn, step 1")
    for iEta in range(eta_rho):
        for iXi in range(xi_rho - 2):
            dx[iEta, iXi + 1] = spheric_dist(
                LON_u[iEta, iXi],
                LAT_u[iEta, iXi],
                LON_u[iEta, iXi + 1],
                LAT_u[iEta, iXi + 1],
            )
    print("dx, dy, pm, pn, step 2")
    for iEta in range(eta_rho):
        dx[iEta, 0] = dx[iEta, 1]
        dx[iEta, xi_rho - 1] = dx[iEta, xi_rho - 2]
    print("dx, dy, pm, pn, step 3")
    for iEta in range(eta_rho - 2):
        for iXi in range(xi_rho):
            dy[iEta + 1, iXi] = spheric_dist(
                LON_v[iEta, iXi],
                LAT_v[iEta, iXi],
                LON_v[iEta + 1, iXi],
                LAT_v[iEta + 1, iXi],
            )
    print("dx, dy, pm, pn, step 4")
    for iXi in range(xi_rho):
        dy[0, iXi] = dy[1, iXi]
        dy[eta_rho - 1, iXi] = dy[eta_rho - 2, iXi]
    print("dx, dy, pm, pn, step 5")
    pm = np.zeros(shape=(eta_rho, xi_rho))
    pn = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            print(
                "iEta=",
                iEta,
                " iXi=",
                iXi,
                " dx=",
                dx[iEta, iXi],
                " dy=",
                dy[iEta, iXi],
            )
            pm[iEta, iXi] = 1 / dx[iEta, iXi]
            pn[iEta, iXi] = 1 / dy[iEta, iXi]
    print("dx, dy, pm, pn, step 6")
    #
    dndx = np.zeros(shape=(eta_rho, xi_rho))
    dmde = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho - 2):
        for iXi in range(xi_rho):
            dmde[iEta + 1, iXi] = 0.5 * (1 / pm[iEta + 2, iXi] - 1 / pm[iEta, iXi])
    for iEta in range(eta_rho):
        for iXi in range(xi_rho - 2):
            dndx[iEta, iXi + 1] = 0.5 * (1 / pn[iEta, iXi + 2] - 1 / pn[iEta, iXi])
    #
    h = np.zeros(shape=(eta_rho, xi_rho))
    f = np.zeros(shape=(eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            f[iEta, iXi] = 2 * (7.29e-5) * sin(LAT_rho[iEta, iXi] * (math.pi / 180))
    #
    mask_rho = np.zeros(shape=(eta_rho, xi_rho))
    angleradMat = get_angle_corr(LON_u, LAT_u)

    GRID_DirectWriteGridFile_Generic(
        GridFile,
        LON_rho,
        LAT_rho,
        LON_u,
        LAT_u,
        LON_v,
        LAT_v,
        LON_psi,
        LAT_psi,
        f,
        angleradMat,
        pm,
        pn,
        dndx,
        dmde,
        xl,
        el,
        mask_rho,
        h,
    )

    LonMin = LON_rho.min()
    LonMax = LON_rho.max()
    LatMin = LAT_rho.min()
    LatMax = LAT_rho.max()
    print("LonMin=", LonMin, " MinLon=", MinLon)
    print("LonMax=", LonMax, " MaxLon=", MaxLon)
    print("LatMin=", LatMin, " MinLat=", MinLat)
    print("LatMax=", LatMax, " MaxLat=", MaxLat)


def WriteNetcdfInitial_Generic(InitFile, eDate, U, V, UBAR, VBAR, ZETA, TEMP, SALT):
    [N, eta_rho, xi_rho] = TEMP.shape
    [eta_u, xi_u] = UBAR.shape
    [eta_v, xi_v] = VBAR.shape
    #
    import netCDF4

    dataset = netCDF4.Dataset(InitFile, "w")
    Nc_s_rho = dataset.createDimension("s_rho", N)
    Nc_s_w = dataset.createDimension("s_w", N + 1)
    Nc_eta_rho = dataset.createDimension("eta_rho", eta_rho)
    Nc_eta_u = dataset.createDimension("eta_u", eta_u)
    Nc_eta_v = dataset.createDimension("eta_v", eta_v)
    Nc_xi_rho = dataset.createDimension("xi_rho", xi_rho)
    Nc_xi_u = dataset.createDimension("xi_u", xi_u)
    Nc_xi_v = dataset.createDimension("xi_v", xi_v)
    Nc_one = dataset.createDimension("one", 1)
    Nc_ocean_time = dataset.createDimension("ocean_time", 0)
    #
    NC_U = dataset.createVariable(
        "u", np.float32, ("ocean_time", "s_rho", "eta_u", "xi_u")
    )
    NC_U.long_name = "u-momentum component"
    NC_U.units = "meter second-1"
    NC_U.time = "ocean_time"
    NC_U.coordinates = "lon_u lat_u s_rho ocean_time"
    NC_U.field = "u-velocity, scalar, series"
    #    import pdb
    #    pdb.set_trace()
    NC_U[0, :, :, :] = U
    #
    NC_V = dataset.createVariable(
        "v", np.float32, ("ocean_time", "s_rho", "eta_v", "xi_v")
    )
    NC_V.long_name = "v-momentum component"
    NC_V.units = "meter second-1"
    NC_V.time = "ocean_time"
    NC_V.coordinates = "lon_v lat_v s_rho ocean_time"
    NC_V.field = "v-velocity, scalar, series"
    NC_V[0, :, :, :] = V
    #
    NC_UBAR = dataset.createVariable(
        "ubar", np.float32, ("ocean_time", "eta_u", "xi_u")
    )
    NC_UBAR.long_name = "vertically integrated u-momentum component"
    NC_UBAR.units = "meter second-1"
    NC_UBAR.time = "ocean_time"
    NC_UBAR.coordinates = "lon_u lat_u ocean_time"
    NC_UBAR.field = "ubar-velocity, scalar, series"
    NC_UBAR[0, :, :] = UBAR
    #
    NC_VBAR = dataset.createVariable(
        "vbar", np.float32, ("ocean_time", "eta_v", "xi_v")
    )
    NC_VBAR.long_name = "vertically integrated v-momentum component"
    NC_VBAR.units = "meter second-1"
    NC_VBAR.time = "ocean_time"
    NC_VBAR.coordinates = "lon_v lat_v ocean_time"
    NC_VBAR.field = "vbar-velocity, scalar, series"
    NC_VBAR[0, :, :] = VBAR
    #
    NC_ZETA = dataset.createVariable(
        "zeta", np.float32, ("ocean_time", "eta_rho", "xi_rho")
    )
    NC_ZETA.long_name = "free-surface"
    NC_ZETA.units = "meter"
    NC_ZETA.time = "ocean_time"
    NC_ZETA.coordinates = "lon_rho lat_rho ocean_time"
    NC_ZETA.field = "free-surface, scalar, series"
    NC_ZETA[0, :, :] = ZETA
    #
    NC_TEMP = dataset.createVariable(
        "temp", np.float32, ("ocean_time", "s_rho", "eta_rho", "xi_rho")
    )
    NC_TEMP.long_name = "potential temperature"
    NC_TEMP.units = "Celsius"
    NC_TEMP.time = "ocean_time"
    NC_TEMP.coordinates = "lon_rho lat_rho s_rho ocean_time"
    NC_TEMP.field = "temperature, scalar, series"
    NC_TEMP[0, :, :] = TEMP
    #
    NC_SALT = dataset.createVariable(
        "salt", np.float32, ("ocean_time", "s_rho", "eta_rho", "xi_rho")
    )
    NC_SALT.long_name = "salinity"
    NC_SALT.units = "PSU"
    NC_SALT.time = "ocean_time"
    NC_SALT.coordinates = "lon_rho lat_rho s_rho ocean_time"
    NC_SALT.field = "salinity, scalar, series"
    NC_SALT[0, :, :] = SALT
    #
    NC_TIME = dataset.createVariable("ocean_time", np.float32, ("ocean_time"))
    NC_TIME.long_name = "time since initialization"
    NC_TIME.units = "seconds since 1968-05-23 00:00:00"
    NC_TIME.calendar = "gregorian"
    NC_TIME.field = "time, scalar, series"
    eTime_start = datetime.datetime(1968, 5, 23, 0, 0)
    delta_time = eDate - eTime_start
    nb_second = delta_time.total_seconds()
    NC_TIME[0] = nb_second
    #
    dataset.close()


def ROMSgetARayVerticalDescription(
    N, Vtransform, Vstretching, Tcline, hc, theta_s, theta_b
):
    Cs_r = N * [0]
    Cs_w = (N + 1) * [0]
    sc_r = N * [0]
    sc_w = (N + 1) * [0]
    half = 1.0 / 2.0
    if Vstretching == 1:
        if theta_s > 0:
            cff1 = 1 / sinh(theta_s)
            cff2 = 1 / (2 * tanh(theta_s / 2))
        else:
            cff1 = -400
            cff2 = -400
        sc_w[0] = -1
        Cs_w[0] = -1
        ds = 1.0 / N
        for k in range(1, N + 1):
            eSc_w = ds * (k - N)
            eSc_r = ds * ((k - N) - half)
            sc_w[k] = eSc_w
            sc_r[k - 1] = eSc_r
            if theta_s > 0:
                Cs_w[k] = (1 - theta_b) * cff1 * sinh(theta_s * eSc_w) + theta_b * (
                    cff2 * tanh(theta_b * (eSc_w + half)) - half
                )
                Cs_r[k - 1] = (1 - theta_b) * cff1 * sinh(theta_s * eSc_r) + theta_b * (
                    cff2 * tanh(theta_b * (eSc_r + half)) - half
                )
            else:
                Cs_w[k] = eSc_w
                Cs_r[k - 1] = eSc_r
    if Vstretching == 2:
        Aweight = 1
        Bweight = 1
        ds = 1.0 / N
        sc_w[N] = 0
        Cs_w[N] = 0
        for k in range(1, N - 1):
            e_sc_w = ds * (k - N)
            sc_w[k] = e_sc_w
            if theta_s > 0:
                Csur = (1 - cosh(theta_s * e_sc_w)) / (cosh(theta_s) - 1)
                if theta_b > 0:
                    Cbot = sinh(theta_s * (e_sc_w + 1)) / sinh(theta_s) - 1
                    Cweight = pow(e_sc_w + 1, Aweight) * (
                        1 + (Aweight / Bweight) * (1 - pow(e_sc_w + 1, Bweight))
                    )
                    Cs_w[k] = Cweight * Csur + (1 - Cweight) * Cbot
                else:
                    Cs_w[k] = Csur
            else:
                Cs_w[k] = e_sc_w
        sc_w[0] = -1
        Cs_w[0] = -1
        for k in range(1, N + 1):
            e_sc_r = ds * ((k - N) - half)
            sc_r[k - 1] = e_sc_r
            if theta_s > 0:
                Csur = (1 - cosh(theta_s * e_sc_r)) / (cosh(theta_s) - 1)
                if theta_b > 0:
                    Cbot = sinh(theta_s * (e_sc_r + 1)) / sinh(theta_s) - 1
                    Cweight = pow(e_sc_r + 1, Aweight) * (
                        1 + (Aweight / Bweight) * (1 - pow(e_sc_r + 1, Bweight))
                    )
                    Cs_r[k - 1] = Cweight * Csur + (1 - Cweight) * Cbot
                else:
                    Cs_r[k - 1] = Csur
            else:
                Cs_r[k - 1] = e_sc_r
    if Vstretching == 3:
        exp_sur = theta_s
        exp_bot = theta_b
        Hscale = 3
        ds = 1.0 / N
        sc_w[N] = 0
        Cs_w[N] = 0
        for k in range(1, N):
            e_sc_w = ds * (k - N)
            sc_w[k] = e_sc_w
            Cbot = log(cosh(Hscale * pow(e_sc_w + 1, exp_bot))) / log(cosh(Hscale)) - 1
            Csur = -log(cosh(Hscale * pow(abs(e_sc_w), exp_sur))) / log(cosh(Hscale))
            Cweight = half * (1 - tanh(Hscale * (e_sc_w + half)))
            Cs_w[k] = Cweight * Cbot + (1 - Cweight) * Csur
        sc_w[0] = -1
        Cs_w[0] = -1
        for k in range(1, N + 1):
            e_sc_r = ds * ((k - N) - half)
            sc_r[k - 1] = e_sc_r
            Cbot = log(cosh(Hscale * pow(e_sc_r + 1, exp_bot))) / log(cosh(Hscale)) - 1
            Csur = -log(cosh(Hscale * pow(abs(e_sc_r), exp_sur))) / log(cosh(Hscale))
            Cweight = half * (1 - tanh(Hscale * (e_sc_r + half)))
            Cs_r[k - 1] = Cweight * Cbot + (1 - Cweight) * Csur
    if Vstretching == 4:
        ds = 1.0 / N
        sc_w[N] = 0
        Cs_w[N] = 0
        for k in range(1, N):
            e_sc_w = ds * (k - N)
            sc_w[k] = e_sc_w
            if theta_s > 0:
                Csur = (1 - cosh(theta_s * e_sc_w)) / (cosh(theta_s) - 1)
            else:
                Csur = -pow(e_sc_w, 2)
            if theta_b > 0:
                Cbot = (exp(theta_b * Csur) - 1) / (1 - exp(-theta_b))
                Cs_w[k] = Cbot
            else:
                Cs_w[k] = Csur
        sc_w[0] = -1
        Cs_w[0] = -1
        for k in range(1, N + 1):
            e_sc_r = ds * ((k - N) - half)
            sc_r[k - 1] = e_sc_r
            if theta_s > 0:
                Csur = (1 - cosh(theta_s * e_sc_r)) / (cosh(theta_s) - 1)
            else:
                Csur = -pow(e_sc_r, 2)
            if theta_b > 0:
                Cbot = (exp(theta_b * Csur) - 1) / (1 - exp(-theta_b))
                Cs_r[k - 1] = Cbot
            else:
                Cs_r[k - 1] = Csur

    return ARVDtyp(
        N, Vstretching, Vtransform, Tcline, hc, theta_s, theta_b, Cs_r, Cs_w, sc_r, sc_w
    )


def GetVertCoord_R(ARVD, hwater, eZeta):
    Vtrans = ARVD.Vtransform
    N = ARVD.N
    z_w = (N + 1) * [0]
    z_r = N * [0]
    Hz = N * [0]
    if Vtrans == 1:
        hinv = 1 / hwater
        z_w[0] = -hwater
        for k in range(1, N + 1):
            cff_r = ARVD.hc * (ARVD.sc_r[k - 1] - ARVD.Cs_r[k - 1])
            cff_w = ARVD.hc * (ARVD.sc_w[k] - ARVD.Cs_w[k])
            cff1_r = ARVD.Cs_r[k - 1]
            cff1_w = ARVD.Cs_w[k]
            z_w0 = cff_w + cff1_w * hwater
            z_w[k] = z_w0 + eZeta * (1 + z_w0 * hinv)
            z_r0 = cff_r + cff1_r * hwater
            z_r[k - 1] = z_r0 + eZeta * (1 + z_r0 * hinv)
            Hz[k - 1] = z_w[k] - z_w[k - 1]
        return VerticalStrat(Hz, z_r, z_w)

    if Vtrans == 2:
        z_w[0] = -hwater
        hinv = 1 / (ARVD.hc + hwater)
        for k in range(1, N + 1):
            cff_r = ARVD.hc * ARVD.sc_r[k - 1]
            cff_w = ARVD.hc * ARVD.sc_w[k]
            cff1_r = ARVD.Cs_r[k - 1]
            cff1_w = ARVD.Cs_w[k]
            cff2_r = (cff_r + cff1_r * hwater) * hinv
            cff2_w = (cff_w + cff1_w * hwater) * hinv
            z_w[k] = eZeta + (eZeta + hwater) * cff2_w
            z_r[k - 1] = eZeta + (eZeta + hwater) * cff2_r
            Hz[k - 1] = z_w[k] - z_w[k - 1]
        return VerticalStrat(Hz, z_r, z_w)

    print("Failed to find relevant entry")
    sys.exit(1)


def CreateStratifiedInitialState(
    InitFile, GridFile, eDate, ARVD, ListDep, ListSalt, ListTemp
):
    import netCDF4

    dataset = netCDF4.Dataset(GridFile, "r")
    DEP = dataset["h"][:, :]
    MSK = dataset["mask_rho"][:, :]
    dataset.close()
    #
    [eta_rho, xi_rho] = DEP.shape
    eta_u = eta_rho
    xi_u = xi_rho - 1
    eta_v = eta_rho - 1
    xi_v = xi_rho
    N = ARVD.N
    s_rho = N
    U = np.zeros(shape=(s_rho, eta_u, xi_u))
    V = np.zeros(shape=(s_rho, eta_v, xi_v))
    UBAR = np.zeros(shape=(eta_u, xi_u))
    VBAR = np.zeros(shape=(eta_v, xi_v))
    ZETA = np.zeros(shape=(eta_rho, xi_rho))
    TEMP = np.zeros(shape=(s_rho, eta_rho, xi_rho))
    SALT = np.zeros(shape=(s_rho, eta_rho, xi_rho))
    hwater = 0
    nb_dep = len(ListDep)

    def get_interpolating_value(ListDep, ListVal, eDep):
        for i in range(1, nb_dep):
            eDep1 = ListDep[i - 1]
            eDep2 = ListDep[i]
            if eDep1 < eDep and eDep <= eDep2:
                alpha2 = (eDep - eDep1) / (eDep2 - eDep1)
                alpha1 = (eDep2 - eDep) / (eDep2 - eDep1)
                eVal1 = ListVal[i - 1]
                eVal2 = ListVal[i]
                return eVal1 * alpha1 + eVal2 * alpha2
        print("Failed to find an entry")
        sys.exit(1)

    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            eDep = DEP[iEta, iXi]
            hwater = eDep
            eZeta = 0
            eVert = GetVertCoord_R(ARVD, hwater, eZeta)
            for iS in range(N):
                eDep = eVert.z_r[iS]
                eTemp = get_interpolating_value(ListDep, ListTemp, eDep)
                eSalt = get_interpolating_value(ListDep, ListSalt, eDep)
                TEMP[iS, iEta, iXi] = eTemp
                SALT[iS, iEta, iXi] = eSalt

    WriteNetcdfInitial_Generic(InitFile, eDate, U, V, UBAR, VBAR, ZETA, TEMP, SALT)


def ComputeDistanceKM(LonDeg1, LatDeg1, LonDeg2, LatDeg2):
    pi = 3.141592653589792
    lon1 = pi * LonDeg1 / (180)
    lat1 = pi * LatDeg1 / (180)
    x1 = cos(lon1) * cos(lat1)
    y1 = sin(lon1) * cos(lat1)
    z1 = sin(lat1)

    lon2 = pi * LonDeg2 / (180)
    lat2 = pi * LatDeg2 / (180)
    x2 = cos(lon2) * cos(lat2)
    y2 = sin(lon2) * cos(lat2)
    z2 = sin(lat2)

    scalprod = x1 * x2 + y1 * y2 + z1 * z2
    EarthRadius = 6371
    if scalprod > 1:
        return 0
    else:
        return EarthRadius * acos(scalprod)
