#include "cfgParser.h"

cfgParser::cfgParser(const string &cfgFileName)
{
    if (fs::exists(cfgFileName))
    {
        cfgFileStream.open(cfgFileName.c_str(), ios::in);
        if (!cfgFileStream.is_open())
        {
            cout << "Cannot open " << cfgFileName << endl;
        }
        string line;
        while (getline(cfgFileStream, line))
        {
            lines.push_back(line);
        }
        parseParameters();
    }
    else
        cout << cfgFileName << " does not exists" << endl;
}

cfgParser::~cfgParser()
{
    if (cfgFileStream.is_open())
        cfgFileStream.close();
}

void cfgParser::splitString(const string &mString)
{
    tokens.clear();
    size_t start = 0;
    size_t end = 0;
    while (end < mString.size())
    {
        if (mString[end] == '\t' || mString[end] == ' ')
        {
            // ignore continuous sep
            if (start < end)
                tokens.push_back(mString.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }
    // save the last token
    tokens.push_back(mString.substr(start));
}

void cfgParser::parseParameters()
{
    for (size_t i = 0; i < lines.size(); i++)
    {
        if (lines[i][0] != '#')
        {
            splitString(lines[i]);
            parametersIXMap.insert({tokens[0], i});
            parametersMap.insert({tokens[0], tokens});
        }
    }
}

size_t cfgParser::findParameter(const string &parameter)
{
    auto iter = parametersIXMap.find(parameter);
    if (iter != parametersIXMap.end())
        return iter->second;
    else
    {
        return 0;
        cout << parameter << "Not Found!" << endl;
    }
}

void cfgParser::setSearch_NameIX()
{
    Search_NameIX = findParameter("Search_Name");
}

void cfgParser::setParent_Mass_WindowsIX()
{
    Parent_Mass_WindowsIX = findParameter("Parent_Mass_Windows");
}

void cfgParser::setElement_PercentIX(const string &element)
{
    string token = "Element_Percent{" + element + "}";
    Element_PercentName = token;
    Element_PercentIX = findParameter(token);
}

void cfgParser::changeSearchName(const string &nameSuffix)
{
    newFileName = parametersMap["Search_Name"][2] + "_" + nameSuffix;
    lines[Search_NameIX] = "Search_Name = " + newFileName;
}

void cfgParser::changeMassWindowsCenter(const int center, const int windowsSize)
{
    vector<int> leftWindows, rightWindows;
    for (int i = 0; i <= windowsSize; i++)
    {
        leftWindows.push_back(center - i);
        rightWindows.push_back(center + i);
    }
    reverse(leftWindows.begin(), leftWindows.end());
    leftWindows.insert(leftWindows.end(), rightWindows.begin() + 1, rightWindows.end());
    string windows;
    for (size_t i = 0; i < leftWindows.size(); i++)
    {
        windows += to_string(leftWindows[i]) + ",";
    }
    windows.pop_back();
    lines[Parent_Mass_WindowsIX] = "Parent_Mass_Windows = " + windows;
}

void cfgParser::changeSIPabundance(const double sipAbundance)
{
    string pct = to_string_with_precision((100.0F - sipAbundance) / 100.0F, 5);
    pct += ",\t" + to_string_with_precision(sipAbundance / 100.0F, 5);
    lines[Element_PercentIX] = Element_PercentName + " \t=\t" + pct;
}

void cfgParser::writeFile(const string &folderPath)
{
    fs::path path{folderPath};
    if (!fs::exists(path))
    {
        cout << path.string() << "Not exists" << endl;
        cout << "Creat" << path.string() << endl;
        fs::create_directories(path);
    }
    path /= newFileName + ".cfg";
    std::ofstream out(path);
    for (string &line : lines)
    {
        out << line << "\n";
    }
    out.close();
}