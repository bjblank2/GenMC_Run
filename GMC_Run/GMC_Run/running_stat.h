#pragma once
class RunningStat
{
public:
    RunningStat();
    void Clear();
    void Push(long double x);
    int NumDataValues();
    long double Mean();
    long double Variance();
private:
    int m_n;
    long double m_oldM, m_newM, m_oldS, m_newS;
};