#include "cfgParser.h"

string configPath;
string outPath;
string element;
string help = "Usage: -i config file path -o config files output path -e element (C or N)";

bool parseArgs(int argc, char const *argv[])
{
    vector<string> vsArguments;
    while (argc--)
        vsArguments.push_back(*argv++);
    for (int i = 1; i < (int)vsArguments.size(); i++)
    {
        if (vsArguments[i] == "-i")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                configPath = vsArguments[i];
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-o")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                outPath = vsArguments[i];
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-e")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                element = vsArguments[i];
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-h" || vsArguments[i] == "--help")
        {
            cout << help << endl;
            return false;
        }
    }
    if (configPath.length() == 0 || outPath.length() == 0 || element.length() == 0)
    {
        cout << help << endl;
        return false;
    }
    return true;
}

bool generateCFGs(string cfgPath, string outPath, string element)
{
    cfgParser parser(cfgPath);
    parser.setSearch_NameIX();
    parser.setParent_Mass_WindowsIX();
    vector<int> centers, widths;
    vector<double> pcts;
    for (size_t pct = 0; pct < 101; pct++)
    {
        pcts.push_back(pct);
    }
    if (element == "C")
    {
        centers = {0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6,
                   7, 7, 8, 9, -2, -1, -1, 0, 0, 1, 1, 2, 3, 3, 4,
                   4, 5, 5, 6, 7, 7, -3, -3, -2, -2, -1, -1, 0, 1,
                   1, 2, 2, 3, 3, 4, 4, 5, -5, -5, -4, -4, -3, -3,
                   -2, -1, -1, 0, 0, 1, 1, 2, 2, 3, 4, -7, -6, -6,
                   -5, -5, -4, -4, -3, -2, -2, -1, -1, 0, 0, 1, 2,
                   2, -8, -8, -7, -7, -6, -6, -5, -4, -4, -3, -3,
                   -2, -2, -1, -1, 0, 0};
        widths = {2, 2, 2, 3, 4, 4, 4, 5, 4, 5, 6, 6, 6, 6, 6,
                  6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 8, 7, 8, 8, 8, 8, 8, 8,
                  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                  8, 8, 8, 8, 8, 7, 8, 7, 7, 7, 7, 6, 7, 6, 6, 6, 6, 6,
                  6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 3, 3, 2};
        pcts[1] = 1.07;
    }
    else if (element == "N")
    {
        centers = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
                   2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5,
                   5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, -1, -5, -5, -5, -5, -4,
                   -4, -4, -4, -4, -4, -4, -4, -4, -3, -3, -3, -3, -3, -3,
                   -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        widths = {3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5,
                  5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 5, 5, 5, 6,
                  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
                  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5,
                  5, 5, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                  5, 5, 5, 4, 4, 4, 4, 3, 3};
        pcts[0] = 0.368;
    }
    else
    {
        cerr << element << " element is not supported!" << endl;
        exit(1);
    }
    parser.setElement_PercentIX(element);
    for (double pct : pcts)
    {
        parser.changeSearchName(to_string_with_precision(pct * 1000.0, 0) + "Pct");
        parser.changeMassWindowsCenter(centers[pct], widths[pct]);
        // set center 0
        // parser.changeMassWindowsCenter(0, widths[pct]);
        parser.changeSIPabundance(pct);
        parser.writeFile(outPath);
    }
    return true;
}

int main(int argc, char const *argv[])
{
    if (parseArgs(argc, argv))
    {
        generateCFGs(configPath, outPath, element);
    }
    return 0;
}
