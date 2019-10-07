#pragma once
class unionFind
{
    //parent relation, parent.put(src,dst) indicates that src points to dst
private:
    int *parents;
    int groups, size;
    bool *isLeader;
public:
    unionFind(int n);

    int find(int src);

    void merge(int v0, int v1);

    int nbOfGroups();

    int *leaders();

    ~unionFind();
};
