#include <iostream>
#include <vector>
#include <cmath>

typedef unsigned int uint;
using namespace std;

template<typename T>
class Field
{
    uint Nx, Ny, Nz;
    T* data_lin;


public:
    uint size(uint dim) const
    {
        if(dim == 0)
            return Nx;
        else if(dim == 1)
            return Ny;
        else if(dim == 2)
            return Nz;
        else
            return 0;
    }


    inline T operator()(uint i, uint j, uint k) const
    {
        return data_lin[k + j * Nz + i * Ny * Nz];
    }
    inline T& operator()(uint i, uint j, uint k)
    {
        return data_lin[k + j * Nz + i * Ny * Nz];
    }

    Field(uint Nx_, uint Ny_, uint Nz_): Nx(Nx_), Ny(Ny_), Nz(Nz_),
        data_lin(new T[Nx*Ny*Nz])
    {
        for(uint i = 0; i < Nx_*Ny_*Nz_; ++i)
            data_lin[i] = 0;
    }

    ~Field()
    {
        delete []data_lin;
    }
};

class WaveguideProfile
{
public:
    virtual bool operator ()(double x, double y, double z) const = 0;
};


class CircularWaveguide : public WaveguideProfile
{
    double radius;
public:
    CircularWaveguide(double radius_):radius(radius_)
    {

    }

    virtual bool operator ()(double x, double y, double z) const
    {
        if(x*x + y*y <= radius*radius)
            return true;
        else
            return false;
    }
};





int main()
{

    const double cl = 299.792458; //mm/ns
    const double pi = 3.1415926;

    double Lx = 6.5; //mm
    double Ly = 6.5; //mm
    double Lz = 60; //mm

    double courant_factor = 0.95;

    double central_frequency = 95; //GHz

    double radius = 2.5;

    double N_lines_transverse = 20;
    double N_lines_Z = 25;


    double central_period = 1 / central_frequency; // ns
    double total_simulation_time = 1000 * central_period;

    //end of the dash panel//

    double central_lambda = cl / central_frequency; // mm



    uint NCellsX = Lx / (central_lambda / N_lines_transverse);
    uint NCellsY = Ly / (central_lambda / N_lines_transverse);
    uint NCellsZ = Lz / (central_lambda / N_lines_Z);

    // Make all the NCells even numbers
    NCellsX = 2 * (NCellsX / 2);
    NCellsY = 2 * (NCellsY / 2);
    NCellsZ = 2 * (NCellsZ / 2);

    cout << "NCells (X, Y, Z):" << NCellsX << ", " << NCellsY << ", " << NCellsZ << endl;

    double dx = Lx / NCellsX;
    double dy = Ly / NCellsY;
    double dz = Lz / NCellsZ;

    double cell_diagonal_length =
            sqrt(dx*dx + dy*dy + dz*dz);

    double dt = courant_factor * cell_diagonal_length / cl; //ns

    cout << "Time period / dt = " << central_period / dt << endl;

    std::vector<double> Xarray(NCellsX + 1);
    std::vector<double> Yarray(NCellsY + 1);
    std::vector<double> Zarray(NCellsZ + 1);

    for(uint i = 0; i <= NCellsX; ++i)
        Xarray[i] = -Lx/2. + Lx * (double)i / (double)NCellsX;

    for(uint i = 0; i <= NCellsY; ++i)
        Yarray[i] = -Ly/2. + Ly * (double)i / (double)NCellsY;

    for(uint i = 0; i <= NCellsZ; ++i)
        Zarray[i] = Lz * (double)i / (double)NCellsZ;


    std::vector<double> XarrayCenters(NCellsX);
    std::vector<double> YarrayCenters(NCellsY);
    std::vector<double> ZarrayCenters(NCellsZ);


    for(uint i = 0; i < NCellsX; ++i)
        XarrayCenters[i] = 0.5 * (Xarray[i] + Xarray[i+1]);

    for(uint i = 0; i < NCellsY; ++i)
        YarrayCenters[i] = 0.5 * (Yarray[i] + Yarray[i+1]);

    for(uint i = 0; i < NCellsZ; ++i)
        ZarrayCenters[i] = 0.5 * (Zarray[i] + Zarray[i+1]);

    Field<bool> MaskCells(NCellsX, NCellsY, NCellsZ);


    WaveguideProfile *profile = new CircularWaveguide(radius);

    for(uint i = 0; i < NCellsX; ++i)
        for(uint j = 0; j < NCellsY; ++j)
            for(uint k = 0; k < NCellsZ; ++k)
            {
                if(profile->operator ()(XarrayCenters[i],
                                        YarrayCenters[j],
                                        ZarrayCenters[k]))
                    MaskCells(i, j, k) = true;
                else
                    MaskCells(i, j, k) = false;
            }

    Field<uint> MaskEx(NCellsX, NCellsY + 1, NCellsZ + 1);
    Field<uint> MaskEy(NCellsX + 1, NCellsY, NCellsZ + 1);
    Field<uint> MaskEz(NCellsX + 1, NCellsY + 1, NCellsZ);


    for(uint i = 0; i < NCellsX; ++i)
        for(uint j = 0; j < NCellsY; ++j)
            for(uint k = 0; k < NCellsZ; ++k)
            {
                if(MaskCells(i, j, k) == false)
                {
                    MaskEx(i, j, k) = 1;
                    MaskEx(i, j+1, k) = 1;
                    MaskEx(i, j, k+1) = 1;
                    MaskEx(i, j+1, k+1) = 1;

                    MaskEy(i, j, k) = 1;
                    MaskEy(i+1, j, k) = 1;
                    MaskEy(i, j, k+1) = 1;
                    MaskEy(i+1, j, k+1) = 1;

                    MaskEz(i, j, k) = 1;
                    MaskEz(i+1, j, k) = 1;
                    MaskEz(i, j+1, k) = 1;
                    MaskEz(i+1, j+1, k) = 1;
                }
            }


    Field<double> Ex(NCellsX, NCellsY + 1, NCellsZ + 1);
    Field<double> Ey(NCellsX + 1, NCellsY, NCellsZ + 1);
    Field<double> Ez(NCellsX + 1, NCellsY + 1, NCellsZ);

    Field<double> Jx(NCellsX, NCellsY + 1, NCellsZ + 1);
    Field<double> Jy(NCellsX + 1, NCellsY, NCellsZ + 1);
    Field<double> Jz(NCellsX + 1, NCellsY + 1, NCellsZ);

    Field<double> Hx(NCellsX + 1, NCellsY, NCellsZ);
    Field<double> Hy(NCellsX, NCellsY + 1, NCellsZ);
    Field<double> Hz(NCellsX, NCellsY, NCellsZ + 1);

    Field<double> Mx(NCellsX + 1, NCellsY, NCellsZ);
    Field<double> My(NCellsX, NCellsY + 1, NCellsZ);
    Field<double> Mz(NCellsX, NCellsY, NCellsZ + 1);

    uint N_time_steps = total_simulation_time / dt;
    double time = 0;

    cout << "Number of time steps = " << N_time_steps << endl;

    for(uint i_time = 0; i_time < N_time_steps; ++i_time)
    {
        if(i_time%10 == 0)
            cout << "Time step " << i_time + 1 << " out of " << N_time_steps << endl;
        //H-field update
        for(uint i = 0; i < NCellsX; ++i)
            for(uint j = 0; j < NCellsY; ++j)
                for(uint k = 0; k < NCellsZ; ++k)
                {
                    Hx(i, j, k) += -dt * cl *
    ( (Ez(i, j+1, k) - Ez(i, j, k)) / dy - (Ey(i, j, k+1) - Ey(i, j, k)) / dz );
                    Hy(i, j, k) += -dt * cl *
    ( (Ex(i, j, k+1) - Ex(i, j, k)) / dz - (Ez(i+1, j, k) - Ez(i, j, k)) / dx );
                    Hz(i, j, j) += -dt * cl *
    ( (Ey(i+1, j, k) - Ey(i, j, k)) / dx - (Ex(i, j+1, k) - Ex(i, j, k)) / dy );

                }



        // Ex update
        for(uint i = 0; i < NCellsX; ++i)
            for(uint j = 1; j < NCellsY - 1; ++j)
                for(uint k = 1; k < NCellsZ - 1; ++k)
                    Ex(i, j, k) += dt * cl *
    ( (Hz(i, j, k) - Hz(i, j-1, k)) / dy - (Hy(i, j, k) - Hy(i, j, k-1)) / dz );


        // Ey update
        for(uint i = 1; i < NCellsX - 1; ++i)
            for(uint j = 0; j < NCellsY; ++j)
                for(uint k = 1; k < NCellsZ - 1; ++k)
                    Ey(i, j, k) += dt * cl *
    ( (Hx(i, j, k) - Hx(i, j, k-1)) / dz - (Hz(i, j, k) - Hz(i-1, j, k)) / dx );

        // Ez update
        for(uint i = 1; i < NCellsX - 1; ++i)
            for(uint j = 1; j < NCellsY - 1; ++j)
                for(uint k = 0; k < NCellsZ; ++k)
                    Ez(i, j, k) += dt * cl *
    ( (Hy(i, j, k) - Hy(i-1, j, k)) / dx - (Hx(i, j, k) - Hx(i, j-1, k)) / dy );


        //Apply mask on E-field

        //Mask on Ex
        for(uint i = 0; i < NCellsX; ++i)
            for(uint j = 0; j < NCellsY + 1; ++j)
                for(uint k = 0; k < NCellsZ + 1; ++k)
                    Ex(i, j, k) *= 1. - MaskEx(i, j, k);

        //Mask on Ey
        for(uint i = 0; i < NCellsX + 1; ++i)
            for(uint j = 0; j < NCellsY; ++j)
                for(uint k = 0; k < NCellsZ + 1; ++k)
                    Ey(i, j, k) *= 1. - MaskEy(i, j, k);

        //Mask on Ez
        for(uint i = 0; i < NCellsX + 1; ++i)
            for(uint j = 0; j < NCellsY + 1; ++j)
                for(uint k = 0; k < NCellsZ; ++k)
                    Ez(i, j, k) *= 1. - MaskEz(i, j, k);

        time += dt;
    }



    cout << "The end" << endl;
    return 0;
}
