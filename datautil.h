
#ifndef PPM_DATAUTIL_H
#define PPM_DATAUTIL_H

#include <cstdio>
#include <memory>
#include <set>
#include "type.h"

using std::set;

namespace PPM::DataUtil
{
std::tuple<vector<int>, Transaction<vector<double>>> read(
    const string file_name, const int event_size);

void binary_label_normalize(vector<int> &alice, vector<int> &bob);
void binary_label_normalize(std::set<int> label_set, vector<int> &label);
void binary_label_normalize(vector<int> &label);

Transaction<Range> discretize(const Transaction<vector<double>> &raw,
                              const vector<vector<double>> &grid,
                              const int event_size);

const vector<vector<double>> partition = {
    {0},
    {-0.4307273, 0.4307273},
    {-0.6744898, 0, 0.6744898},
    {-0.841621233572914, -0.2533471031358, 0.2533471031358, 0.841621233572914},
    {-0.967421566101701,
     -0.430727299295457,
     0,
     0.430727299295457,
     0.967421566101701},
    {-1.06757052387814,
     -0.565948821932863,
     -0.180012369792705,
     0.180012369792705,
     0.565948821932863,
     1.06757052387814},
    {-1.15034938037601,
     -0.674489750196082,
     -0.318639363964375,
     0,
     0.318639363964375,
     0.674489750196082,
     1.15034938037601},
    {-1.22064034884735,
     -0.764709673786387,
     -0.430727299295457,
     -0.139710298881862,
     0.139710298881862,
     0.430727299295457,
     0.764709673786387,
     1.22064034884735},
    {-1.2815515655446,
     -0.841621233572914,
     -0.524400512708041,
     -0.2533471031358,
     0,
     0.2533471031358,
     0.524400512708041,
     0.841621233572914,
     1.2815515655446},
    {-1.33517773611894,
     -0.908457868537385,
     -0.604585346583237,
     -0.348755695517045,
     -0.114185294321428,
     0.114185294321428,
     0.348755695517045,
     0.604585346583237,
     0.908457868537385,
     1.33517773611894},
    {-1.38299412710064,
     -0.967421566101701,
     -0.674489750196082,
     -0.430727299295457,
     -0.210428394247925,
     0,
     0.210428394247925,
     0.430727299295457,
     0.674489750196082,
     0.967421566101701,
     1.38299412710064},
    {-1.42607687227285,
     -1.0200762327862,
     -0.736315917376129,
     -0.502402223373355,
     -0.293381232121193,
     -0.0965586152896391,
     0.0965586152896394,
     0.293381232121194,
     0.502402223373355,
     0.73631591737613,
     1.0200762327862,
     1.42607687227285},
    {-1.46523379268552,
     -1.06757052387814,
     -0.791638607743375,
     -0.565948821932863,
     -0.36610635680057,
     -0.180012369792705,
     0,
     0.180012369792705,
     0.36610635680057,
     0.565948821932863,
     0.791638607743375,
     1.06757052387814,
     1.46523379268552},
    {-1.50108594604402,
     -1.11077161663679,
     -0.841621233572914,
     -0.622925723210088,
     -0.430727299295457,
     -0.2533471031358,
     -0.0836517339071291,
     0.0836517339071291,
     0.2533471031358,
     0.430727299295457,
     0.622925723210088,
     0.841621233572914,
     1.11077161663679,
     1.50108594604402},
    {-1.53412054435255,
     -1.15034938037601,
     -0.887146559018876,
     -0.674489750196082,
     -0.488776411114669,
     -0.318639363964375,
     -0.157310684610171,
     0,
     0.157310684610171,
     0.318639363964375,
     0.488776411114669,
     0.674489750196082,
     0.887146559018876,
     1.15034938037601,
     1.53412054435255},
    {-1.5647264713618,
     -1.18683143275582,
     -0.928899491647271,
     -0.721522283982343,
     -0.541395085129088,
     -0.377391943828554,
     -0.223007830940367,
     -0.0737912738082727,
     0.0737912738082727,
     0.223007830940367,
     0.377391943828554,
     0.541395085129088,
     0.721522283982343,
     0.928899491647271,
     1.18683143275582,
     1.5647264713618},
    {-1.59321881802305,
     -1.22064034884735,
     -0.967421566101701,
     -0.764709673786387,
     -0.589455797849779,
     -0.430727299295457,
     -0.282216147062508,
     -0.139710298881862,
     0,
     0.139710298881862,
     0.282216147062508,
     0.430727299295457,
     0.589455797849779,
     0.764709673786387,
     0.967421566101701,
     1.22064034884735,
     1.59321881802305},
    {-1.61985625863827,
     -1.25211952026522,
     -1.00314796766253,
     -0.8045963803603,
     -0.633640000779701,
     -0.47950565333095,
     -0.336038140371823,
     -0.199201324789267,
     -0.0660118123758407,
     0.0660118123758406,
     0.199201324789267,
     0.336038140371823,
     0.47950565333095,
     0.633640000779701,
     0.8045963803603,
     1.00314796766253,
     1.25211952026522,
     1.61985625863827},
    {-1.64485362695147,
     -1.2815515655446,
     -1.03643338949379,
     -0.841621233572914,
     -0.674489750196082,
     -0.524400512708041,
     -0.385320466407568,
     -0.2533471031358,
     -0.125661346855074,
     0,
     0.125661346855074,
     0.2533471031358,
     0.385320466407568,
     0.524400512708041,
     0.674489750196082,
     0.841621233572914,
     1.03643338949379,
     1.2815515655446,
     1.64485362695147},
    {-1.6683911939470795,  -1.309171716785777,   -1.0675705238781414,
     -0.8761428492468408,  -0.712443032389489,   -0.5659488219328631,
     -0.43072729929545756, -0.30298044805620655, -0.1800123697927051,
     -0.05971709978532289, 0.05971709978532289,  0.18001236979270496,
     0.30298044805620655,  0.43072729929545744,  0.5659488219328631,
     0.7124430323894889,   0.8761428492468408,   1.0675705238781412,
     1.309171716785777,    1.668391193947079},
    {// 22
     -1.6906216295848977,
     -1.335177736118937,
     -1.0968035620935135,
     -0.9084578685373851,
     -0.7478585947633022,
     -0.6045853465832371,
     -0.4727891209922674,
     -0.3487556955170447,
     -0.22988411757923205,
     -0.11418529432142838,
     0.0,
     0.11418529432142824,
     0.2298841175792322,
     0.3487556955170447,
     0.4727891209922672,
     0.6045853465832371,
     0.7478585947633022,
     0.9084578685373853,
     1.0968035620935135,
     1.3351777361189363,
     1.6906216295848986},
    {// 23
     -1.7116753065097285, -1.3597373839386062,  -1.1243382315686392,
     -0.9388143168769032, -0.7810338115227088,  -0.6406668899191049,
     -0.5119362138713294, -0.39119625818947174, -0.2759210631079599,
     -0.164210777079331,  -0.05451891484810106, 0.05451891484810106,
     0.16421077707933085, 0.27592106310796005,  0.39119625818947174,
     0.5119362138713294,  0.6406668899191048,   0.7810338115227088,
     0.9388143168769032,  1.1243382315686385,   1.3597373839386055,
     1.7116753065097288},
    {// 24
     -1.731664396122245,  -1.382994127100638,   -1.1503493803760079,
     -0.967421566101701,  -0.8122178014999129,  -0.6744897501960817,
     -0.5485222826980979, -0.43072729929545756, -0.31863936396437514,
     -0.2104283942479247, -0.1046334556140754,  0.0,
     0.10463345561407526, 0.21042839424792484,  0.31863936396437514,
     0.43072729929545744, 0.5485222826980981,   0.6744897501960817,
     0.8122178014999129,  0.967421566101701,    1.1503493803760079,
     1.382994127100638,   1.7316643961222453},
    {// 25
     -1.75068607125217,   -1.4050715603096329,  -1.1749867920660904,
     -0.994457883209753,  -0.8416212335729142,  -0.7063025628400874,
     -0.5828415072712162, -0.46769879911450823, -0.3584587932511938,
     -0.2533471031357997, -0.15096921549677725, -0.05015358346473367,
     0.05015358346473367, 0.1509692154967774,   0.2533471031357997,
     0.3584587932511938,  0.4676987991145084,   0.5828415072712162,
     0.7063025628400874,  0.8416212335729143,   0.994457883209753,
     1.1749867920660904,  1.4050715603096329,   1.7506860712521692}};

}  // namespace PPM::DataUtil

#endif
