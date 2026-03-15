bool cEOP::ReadEOP(string directory) {
    // 自动处理路径结尾的斜杠
    if (directory.back() != '/' && directory.back() != '\\') directory += "/";
    
    string eopFile = directory + "EOP-All.txt";
    string taiFile = directory + "tai-utc.dat";

    ifstream inEOP(eopFile);
    if (!inEOP.is_open()) return false;
    
    m_pstEOP.clear();
    m_nNumberEOP = 0;
    string line;
    while (getline(inEOP, line)) {
        if (line.empty() || line[0] == '#') continue;
        EOPData eop;
        stringstream ss(line);
        ss >> eop.dfMJD >> eop.dfX >> eop.dfY >> eop.dfUT1UTC;
        // 假设简化读取，只读前4个
        eop.dfJD = eop.dfMJD + 2400000.5; // 近似
        eop.dfX *= 4.84813681109536e-6; // arcsec to rad
        eop.dfY *= 4.84813681109536e-6;
        m_pstEOP.push_back(eop);
        m_nNumberEOP++;
    }
    inEOP.close();

    ifstream inTAI(taiFile);
    if (!inTAI.is_open()) return false;
    
    m_nNumberTAIUTC = 0;
    while (getline(inTAI, line)) {
        if (line.empty() || line[0] == '#') continue;
        stTAIUTC tai;
        stringstream ss(line);
        ss >> tai.dfJDUTC >> tai.dfJDTAI >> tai.df0 >> tai.dfMJD >> tai.df1;
        m_pstTAIUTC[m_nNumberTAIUTC] = tai;
        m_nNumberTAIUTC++;
        if (m_nNumberTAIUTC >= 1000) break;
    }
    inTAI.close();

    m_bInit = true;
    return true;
}