#pragma once
class cell
{
private:
    int leader;
    double min;

public:
    cell(int leader_);
    int getLeader();
    double getMin();
    void setMin(double min);
    ~cell();
};