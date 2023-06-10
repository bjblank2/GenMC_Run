#pragma once
class RunningStat
{
public:
    RunningStat();
    void Clear();
    void Push(float x);
    int NumDataValues();
    float Mean();
    float Variance();
private:
    int m_n;
    float m_oldM, m_newM, m_oldS, m_newS;
};