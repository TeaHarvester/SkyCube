#include "Map.h"

void Map::Write()
{
    std::ofstream file1;
    std::ofstream file2;
    file1.open("topology.txt");
    file2.open("glacier.txt");
    std::vector<int>::iterator iter1;
    std::vector<int>::iterator iter2 = snowpack.begin();;

    for (iter1 = terrain.begin(); iter1 < terrain.end(); ++iter1, ++iter2)
    {
        file1 << *iter1 << std::endl;
        file2 << *iter2 << std::endl;
    }

    file1.close();
    file2.close();
}

void Map::Faultline()
{
    std::cout << std::endl << "generating faultline...";

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
        std::uniform_int_distribution<> height_dist(faultlines[(*iter)], faultlines[(*iter)]*2);
        std::uniform_int_distribution<> base_dist(faultlines[(*iter)], faultlines[(*iter)]*2);

        unsigned int height = (int)std::round(height_dist(gen));
        unsigned int base = (int)std::round(base_dist(gen));
        Mountain(*iter, height, base);   
    }
}

void Map::DrawTerrain()
{
    float mean_height = (float)std::accumulate(terrain.begin(), terrain.end(), 0)/(float)terrain.size();
    std::vector<int> iter;

    for (unsigned int i = 0; i < terrain.size(); ++i)
    {
        if ((float)terrain[i] > 1.5*mean_height)
        {
            features[i] = '^';
        }

        else if (runoff[i])
        {
            features[i] = '~';
        }

        else if ((float)terrain[i] > 0.8*mean_height)
        {
            features[i] = '.';
        }

        else 
        {
            features[i] = ' ';
        }

        // if (snowpack[i] && terrain[i] > temperature)
        // {
        //     features[i] = '*';
        // }
    }
}

std::vector<int> Map::Cut(const int centre, bool river, std::vector<int>& upstream)
{
    std::vector<int> output;
    const int xcentre = centre/ysize;
    const int ycentre = centre%ysize;
    const int st = 2*stagger[ycentre] - 1;
    int direction[6];
    direction[0] = centre - ysize;
    direction[1] = centre + ysize;
    direction[2] = centre - 1;
    direction[3] = centre - 1 + st*ysize;
    direction[4] = centre + 1;
    direction[5] = centre + 1 + st*ysize;

    std::vector<int> grad(3, 0);
    std::vector<int> flow(3, -10);

    if (xcentre - 1 >= 0 && xcentre + 1 < xsize)
    {
        grad[0] = terrain[direction[1]] - terrain[direction[0]];
        
        if (grad[0] > 0)
        {
            flow[0] = direction[0];
        }

        else if (grad[0] < 0)
        {
            flow[0] = direction[1];
        }

        if (ycentre - 1 >= 0 && ycentre + 1 < ysize)
        {
            grad[1] = terrain[direction[4]] - terrain[direction[3]];
            grad[2] = terrain[direction[5]] - terrain[direction[2]];

            if (grad[1] > 0)
            {
                flow[1] = direction[3];
            }

            else if (grad[1] < 0)
            {
                flow[1] = direction[4];
            }

            if (grad[2] > 0)
            {
                flow[2] = direction[2];
            }

            else if (grad[2] < 0)
            {
                flow[2] = direction[5];
            }
        }

        else
        {
            grad.erase(grad.end() - 2, grad.end());
            flow.erase(flow.end() - 2, flow.end());
        }

        std::vector<int>::iterator origin;

        for (origin = upstream.begin(); origin < upstream.end(); ++origin)
        {
            std::vector<int>::iterator grad_iter = grad.begin();
            std::vector<int>::iterator flow_iter;

            for (flow_iter = flow.begin(); flow_iter < flow.end(); ++flow_iter, ++grad_iter)
            {
                if (*flow_iter == *origin)
                {
                    flow.erase(flow_iter, flow_iter + 1);
                    grad.erase(grad_iter, grad_iter + 1);
                    break;
                }
            }
        }

        int max_grad = 0;
        std::vector<int>::iterator grad_iter;

        for (grad_iter = grad.begin(); grad_iter < grad.end(); ++grad_iter)
        {
            int grad_value = std::abs(*grad_iter);

            if (grad_value > max_grad)
            {
                max_grad = grad_value;
            }
        }

        grad_iter = grad.begin();
        std::vector<int>::iterator flow_iter;

        for (flow_iter = flow.begin(); flow_iter < flow.end(); ++flow_iter, ++grad_iter)
        {
            if (std::abs(*grad_iter) == max_grad)
            {
                output.push_back(*flow_iter);
            }
        }

        if (output.size() == 3)
        {
            std::random_device rd; 
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> distrib(0, 1);
            int rotation = distrib(gen);

            if (rotation)
            {
                output[0] = direction[0];                
                output[1] = st == 1 ? direction[3] : direction[2];
                output[2] = st == 1 ? direction[5] : direction[4];
            }

            else
            {
                output[0] = direction[1];                
                output[1] = st == 1 ? direction[2] : direction[3];
                output[2] = st == 1 ? direction[4] : direction[5];
            }
        }
    }

    std::vector<int>::iterator cutter;

    for (cutter = output.begin(); cutter < output.end(); ++cutter)
    {
        if (river)
        {
            ++runoff[*cutter];
        }

        else if (terrain[*cutter] > 0)
        {
            ++snowpack[*cutter];
            --terrain[*cutter];

            // if (terrain[*cutter] > 0)
            // {
            //     --terrain[*cutter];
            // }
        }
    }

    cutter = output.begin();

    while (cutter < output.end())
    {
        if (!terrain[*cutter] && !river)
        {
            output.erase(cutter, cutter + 1);
            cutter = output.begin();
            continue;
        }

        ++cutter;
    }

    return output;
}

void Map::Glaciate(const int height, bool river)
{
    std::cout << std::endl << "glaciating...";

    if (sources.begin() == sources.end())
    {
        for (int i = 0; i < xsize*ysize; ++i)
        {
            if (terrain[i] >= height)
            {
                sources.push_back(i);
                ++snowpack[i];
            }
        }
    }

    std::vector<int>::iterator source;

    for (source = sources.begin(); source < sources.end(); ++source)
    {
        std::vector<int> cuttable(source, source + 1);
        std::vector<int>::iterator cutter;
        std::vector<int> immune(1, -1);

        while (cuttable.size() > 0)
        {
            cutter = cuttable.end() - 1;
            std::vector<int> cutted = Cut(*cutter, river, immune);
            cuttable.pop_back();
            cuttable.insert(cuttable.end(), cutted.begin(), cutted.end());
            immune.insert(immune.end(), cutted.begin(), cutted.end());
        }
    }
}

void Map::Mountain(const int centre, const float peak, const float width)
{
    terrain[centre] = terrain[centre] < peak ? peak : terrain[centre];
    const int xcentre = centre/ysize;
    const int ycentre = centre%ysize;
    int base = ceil((width + 1)/2);
    const float gradient = 2*peak/width;
    std::vector<int> elevated;
    std::vector<float> elevation;

    for (int i = 0; i < base; ++i)
    {
        if (xcentre - (i + 1) >= 0)
        {
            elevated.push_back(ycentre + (xcentre - (i + 1))*ysize);
            elevation.push_back(std::round(peak - i*gradient));
        }

        else
        {
            int xwrapper = xsize + xcentre - (i + 1);
            elevated.push_back(ycentre + xwrapper*ysize);
            elevation.push_back(std::round(peak - (i + 1)*gradient));
        }

        elevated.push_back(ycentre + ((xcentre + (i + 1))%xsize)*ysize);
        elevation.push_back(std::round(peak - (i + 1)*gradient));

        int ntile = width - 1;
        int ycoord1 = ycentre - (i + 1);
        int ycoord2 = ycentre + (i + 1);

        for (int j = 0; j < ntile; ++j)
        {
            int xcoord = xcentre - ntile/2 + stagger[ycentre] + j;

            if (xcoord >= xsize)
            {
                xcoord -= xsize;
            }

            else if (xcoord < 0)
            {
                xcoord += xsize;
            }

            if (ycoord2 >= ysize)
            {
                ycoord2 -= ysize;
            }

            if (ycoord1 < 0)
            {
                ycoord1 += ysize;
            }



            float xdist = (float)(j - ntile/2) + (float)stagger[ycoord1]/2;
            float ydist = sqrt(pow((float)(i + 1), 2) - pow(0.5, 2));
            float dist = sqrt(pow(xdist, 2) + pow(ydist, 2));

            if (ycoord1 >= 0)
            {
                elevated.push_back(ycoord1 + xcoord*ysize);
                elevation.push_back(std::round(peak - dist*gradient));
            }

            else 
            {
                int ywrapper = ysize + ycoord1;
                elevated.push_back(ywrapper + xcoord*ysize);
                elevation.push_back(std::round(peak - dist*gradient));
            }

            elevated.push_back(ycoord2%ysize + xcoord*ysize);
            elevation.push_back(std::round(peak - dist*gradient));
        }

        for (unsigned int i = 0; i < elevated.size(); ++i)
        {
            if (terrain[elevated[i]] < elevation[i])
            {
                terrain[elevated[i]] = elevation[i];
            }
        }
    }
}

std::vector<int> Map::GetDir(const int centre)
{
    const int xcentre = centre/ysize;
    const int ycentre = centre%ysize;

    std::vector<int> dir(4);
    dir[0] = xcentre > 0 ? ycentre + ysize*((xcentre - 1)) : ycentre + ysize*(xsize - 1);
    dir[1] = xcentre < xsize - 1 ? ycentre + ysize*((xcentre + 1)) : ycentre;
    dir[2] = ycentre > 0 ? centre - 1 : xcentre*ysize + (ysize - 1);
    dir[3] = ycentre < ysize - 1 ? centre + 1 : xcentre*ysize;

    return dir;
}

void Map::Print(bool topology)
{
    std::cout << std::endl;

    for (int i = 0; i < ysize; ++i)
    {
        std::cout << std::endl;

        if (stagger[i])
        {
            std::cout << ' ';
        }

        for (int j = 0; j < xsize; ++j)
        {
            if (topology)
            {
                std::cout << terrain[i + j*ysize] << ' ';
            }

            else 
            {
                std::cout << faultlines[i + j*ysize] << ' ';
            }
        }
    }
}

Map::Map(int sx, int sy)
:
terrain(sx*sy, 0),
atmos(NULL),
xsize(sx),
ysize(sy),
runoff(sx*sy, 0),
snowpack(sx*sy, 0),
sources(),
faultlines(sx*sy, 0),
features(sx*sy, ' '),
stagger(sy, 0)
{
    for (int i = 0; i < ysize; i = i + 2)
    {
        stagger[i] = 1;
    }

    Faultline();
    // Glaciate(16, false);
    // terrain[size*(size/2) + (size/2)] += 1;
    // Mountain(size*(size/2) + (size/2), 3, 5);
    // Write();
    Print(true);
    atmos = new Atmosphere(xsize, ysize, 6, 10000, 101500, 50000, 60, terrain);
}

Map::Map(const int nlat)
:
terrain(6*nlat*nlat, 0),
xsize(nlat),
ysize(nlat),
runoff(6*nlat*nlat, 0),
snowpack(6*nlat*nlat, 0),
sources(),
faultlines(6*nlat*nlat, 0),
features(6*nlat*nlat, ' '),
stagger(6*nlat*nlat, 0)
{
}

Map::~Map()
{
    if (atmos)
    {
        delete atmos;
        atmos = NULL;
    }
}