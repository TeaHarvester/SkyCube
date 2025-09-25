#include "Earth.h"

void Earth::Faultline(int tile)
{
    std::cout << std::endl << "generating faultline...";


    int t_ = tile*nlat*nlat; 
    std::random_device rd; 
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> xdistrib(ceil((float)xsize/3), floor(2*(float)xsize/3));
    std::uniform_int_distribution<> ydistrib(ceil((float)ysize/3), floor(2*(float)ysize/3));
    std::uniform_int_distribution<> markov(0, 3);

    const int true_xcentre = xdistrib(gen);
    const int true_ycentre = ydistrib(gen);

    int xcentre = true_xcentre;
    int ycentre = true_ycentre;
    int centre = ysize*xcentre + ycentre;

    std::vector<int> direction = GetDir(centre);
    
    bool first_pass = true;
    bool map_bounds = xcentre > 0 && xcentre < xsize - 1 &&
                      ycentre > 0 && ycentre < ysize - 1;

    std::vector<int> fault(1, centre);

    while (map_bounds)
    {
        int idir = markov(gen);
        int destination = direction[idir];

        std::vector<int>::iterator iter;
        int breakout = 0;

        for (iter = fault.begin(); iter < fault.end(); ++iter)
        {
            ++breakout;

            if (breakout > 20)
            {
                break;
            }

            if (destination == *iter)
            { 
                iter = fault.begin();
                idir = markov(gen);
                destination = direction[idir];
                fault.pop_back();
                iter = fault.begin();
                continue;
            }
        }

        centre = destination;
        xcentre = centre/ysize;
        ycentre = centre%ysize;
        direction = GetDir(centre);
        fault.push_back(centre);
        faultlines[centre] += 1;

        map_bounds = xcentre > 0 && xcentre < xsize - 1 &&
                     ycentre > 0 && ycentre < ysize - 1;

        if (!map_bounds && first_pass)
        {
            first_pass = false;
            xcentre = true_xcentre;
            ycentre = true_ycentre;
            centre = ysize*xcentre + ycentre;
            direction = GetDir(centre);
            map_bounds = true;
        }
    }


    std::vector<int>::iterator iter;

    for (iter = fault.begin(); iter < fault.end(); ++iter)
    {
        int activity = faultlines[(*iter)];
        std::uniform_int_distribution<> height_dist(activity/2, 2*activity/2);
        std::uniform_int_distribution<> base_dist(activity/2, 2*activity/2);

        unsigned int height = (int)std::round(height_dist(gen));
        unsigned int base = (int)std::round(base_dist(gen));
        Mountain(*iter + t_, height, base);   
    }

    faultlines = std::vector<int>(6*nlat*nlat, 0);
}

void Earth::Print(int tile)
{
    int tsize = nlat*nlat;
    int t_ = tile*tsize;

    std::cout << std::endl;

    for (int i = 0; i < ysize; ++i)
    {
        std::cout << std::endl;

        for (int j = 0; j < xsize; ++j)
        {
            std::cout << terrain[i + j*ysize + t_] << ' ';
        }
    }
}

Earth::Earth(const int nlat_) : Map(nlat_),
nlat(nlat_)
{
    for (int i = 0; i < 6; ++i)
    {
        Faultline(i);
        Print(i);
    }
}