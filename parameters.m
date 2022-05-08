%% Class parameters.m
% This class contains the more relevant physical parameters of the project.

classdef parameters
    properties (Constant)
        perihelion =  2.06650E+08; % [km] Perihelion
        aphelion =    2.49261E+08; % [km] Aphelion 
        a =           2.28E+08; % [km] Semi-major axis
        T =           686.973*3600; % [h] Orbit period
        Ts =          779.94*3600; % [h] Synodic period
        i =           1.848; % [°] Orbit inclination
        e =           0.0935; % [-] Eccentricity
        SidRot =      24.6229*3600; % [s] Sideral rotation period
        Tday =        24.6597*3600; % [s] Length of day
        M =           6.4169E+23; % [kg] Mass
        V =           1.6318E+11; % [km^3] Volume
        R =           3396.2; % [km] Equatorial radius
        Rpolar =      3376.2; % [km] Polar radius
        VMR =         3389.5; % [km] Volumetric Mean Radius
        Vcore =       1700; % [km] Core radius
        ellipticity = 0.00589; % [km] Ellipticity (flattening)
        rhoMean =     3934; % [kg/m^3] Mean density
        gSurf =       3.71; % [m/s^2] Surface gravity
        aSurf =       3.69; % [m/s^2] Surface acceleration
        vEscape =     5.03; % [km/s] Escape velocity
        mu =          4.2828E+04; % [km^3/s^2] Viscosity
        albedo =      0.15; % [-] Albedo
        vband =       -1.6; % [V] V-band magnitude
        Js =          590.2792; % [W/m^2] Solar irradiance
        Tb =          209.8; % [K] Black-body temperature
        p =           6.36; % [Mbar] Surface pressure
        rhoSurf =     0.02; % [kg/m^3] Surface density
        Tmean =       210; % [K] Average T
        D =           228e6; % [km] Distance Sun-Mars [https://solarsystem.nasa.gov/planets/mars]
        eps =         1; % [-] Emissivity of Mars
        P =           3.856e26; % [W] Total power output from the Sun [Spacecraft System Engineering]
        
        % Physical constants
        sigma =       5.67e-8; % [Wm^-2K^-4] Stefan-Boltzmann constant
    end
end