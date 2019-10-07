#include "cell.hpp"

cell::cell(int leader_){
    leader = leader_;
    min = -1;
}
int cell::getLeader(){
    return leader;
}
double cell::getMin(){
    return min;
}
void cell::setMin(double min_){
    min = min_;
}
cell::~cell(){

}