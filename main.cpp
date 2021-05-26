#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <set>
#include <ctime>
#include <algorithm>
#include <chrono>

struct coordStruct {
    int x;
    int y;
};

struct part {
    int vert;
    double weight;
};

struct distStruct {
    double g;
    double rhs;
};

struct forQueue {
    int v;
    double key1;
    double key2;
};

struct Compare
{
    bool operator()(const forQueue& x, const forQueue& y) const
    {
        if (x.key1 != y.key1) {
            return x.key1 > y.key1;
        } else {
            return x.key2 > y.key2;
        }
    }
};

using namespace std;
using namespace std::chrono;

int vertNum, edgeNum, startVert, endVert;

priority_queue<forQueue, vector<forQueue>, Compare> Q;
vector <vector <part>> list(0);
vector <distStruct> dist (0);
vector <coordStruct> coord(0);
int fromVert = 0;
const auto MAX = 100000;
const auto eps = 0.001;
int counter = 0;
high_resolution_clock::time_point start_time;
high_resolution_clock::time_point end_time;
duration<double> fullDijTime;


double Dijkstra() {
    //startDijTime = clock();
    //startDijTime = high_resolution_clock::now();
    vector <double> distDijkstra(vertNum+1);
    for (int i = 0; i <= vertNum; i++) {
        distDijkstra[i] = INT_MAX;
    }
    distDijkstra[fromVert] = 0;
    set <pair <int, double> > s;
    s.insert(pair<int, double>(fromVert, 0));
    while (!s.empty()) {
        int from = s.begin()->first;
        s.erase(s.begin());

        for (int i = 0; i < list[from].size(); i++) {
            part to = list[from][i];

            if (distDijkstra[from] + to.weight < distDijkstra[to.vert]) {
                s.erase(pair<int, double>(to.vert, to.weight));
                distDijkstra[to.vert] = distDijkstra[from] + to.weight;
                s.insert(pair<int, double>(to.vert, distDijkstra[to.vert]));
            }
        }
    }
    //endDijTime = clock();
    //endDijTime = high_resolution_clock::now();
    //fullDijTime += endDijTime - startDijTime;
    //fullDijTime += duration_cast<duration<double>>(endDijTime - startDijTime);
    if (distDijkstra[endVert] == INT_MAX) {
        return(-1);
    }
    else {
        return(distDijkstra[endVert]);
    }
}

void clearData() {
    while (!Q.empty()) {
        Q.pop();
    }
    list.clear();
    list.resize(20);
    dist.clear();
    coord.clear();
}

bool inVerts(vector<int> vec, int s, int e) {
    for (auto i: vec) {
        if (i == s || i == e) {
            return true;
        }
    }
    return false;
}

int operator< (forQueue a, forQueue b) { //tested
    if (a.key1 < b.key1) {
        return true;
    }
    if (a.key1 > b.key1) {
        return false;
    }
    if (a.key2 <= b.key2) {
        return true;
    }
    return false;
}

vector<int> pred(int u) { //tested
    vector<int> res;
    for (int i = 0; i < list.size(); i++) {
        for (auto j: list[i]) {
            if (j.vert == u) {
                res.push_back(i);
            }
        }
    }
    return res;
}

vector<int> succ(int f) { //tested
    vector<int> res;
    for (auto j: list[f]) {
        res.push_back(j.vert);
    }
    return res;
}

double getEdgeValue(int a, int b) { //tested
    for (auto j: list[a]) {
        if (j.vert == b) {
            return j.weight;
        }
    }
}

void remove(int u) { //tested
    vector<forQueue> vec;
    while (!Q.empty()) {
        if (Q.top().v != u) {
            vec.push_back(Q.top());
        }
        Q.pop();
    }
    for (auto i: vec) {
        Q.push(i);
    }
}

double h(int a, int b) { //евклидово расстояние по прямой
    return sqrt(pow(coord[a].x - coord[b].x, 2) + pow(coord[a].y - coord[b].y, 2));
}

double minSuccVert(int prevVert, char val) { //return path or vert
    double min = MAX;
    double pathToS = 0;
    int vert = 0;
    for (auto s : succ(prevVert)) {
        pathToS = getEdgeValue(prevVert, s) + dist[s].g;
        if (min > pathToS) {
            min = pathToS;
            vert = s;
        }
    }
    if (val == 'v') {
        return vert;
    } else {
        return min;
    }
}

bool inQueue(int u) { //tested
    priority_queue<forQueue, vector<forQueue>, Compare> copyQ;
    copyQ = Q;
    while (!copyQ.empty()) {
        if (copyQ.top().v == u) {
            return true;
        }
        copyQ.pop();
    }
    return false;
}

forQueue calcKey(int s) {
    forQueue mem{};
    mem.v = s;
    mem.key1 = min(dist[s].g, dist[s].rhs) + h(startVert, s);
    mem.key2 = min(dist[s].g, dist[s].rhs);
    return mem;
}

void initialize() {
    for (auto i = 0; i < vertNum; i++) {
        dist.push_back({MAX, MAX});
    }
    dist[endVert].rhs = 0;
    Q.push(calcKey(endVert));
}

void updateVertex(int u) {
    if (u != endVert) {
        dist[u].rhs = minSuccVert(u, 'p');
    }

    if (inQueue(u)) {
        remove(u);
    }

    if (dist[u].g != dist[u].rhs) {
        Q.push(calcKey(u));
    }
}

void computeShortestPath() {
    while (((Q.top() < calcKey(fromVert)) ||
        (dist[startVert].rhs != dist[startVert].g)) && (!Q.empty()))   {
        auto u = Q.top().v;
        Q.pop();
        if (dist[u].g > dist[u].rhs) {
            dist[u].g = dist[u].rhs;
            for (auto s : pred(u)) {
                updateVertex(s);
            }
        } else {
            dist[u].g = MAX;
            for (auto s : pred(u)) {
                updateVertex(s);
            }
            updateVertex(u);                                                                                                                           //may cause troubles
        }
    }
}

int main() {
    //input : vert number, edge number, from, to
    //vert number lines with x, y coordinates
    //edge number lines with a vert, b vert and edge value

    //if changed
    //Y
    //changed vert number, changed edge number
    //changed vert number with vert, x, y
    //changed edge number with start, end, changed edge value

    //0-10 - 5vert
    //11-20 - 10vert
    //21-25 - 100vert
    //26-40 - 1000vert
    //41-60 - 5000vert
    //60-64 - 5vert with change
    //65-79 - 10vert with change
    //80-99 - 100vert with change
    //100-109 - 1000vert with change
    //110-119 - 5000vert with change
    int testNum = 0;
    double timeD = 0;
    vector<duration<double>> times;
    for (int i = 1; i <= 3; i++) {
    //for (int j = 0; j < 10; j++) {
        //9 тест быстрее

            clearData();
            testNum = i;
            string testName = "input" + to_string(testNum) + ".txt";
            string output = "output" + to_string(testNum) + ".txt";
            ifstream in("D:\\6sem\\algorithm\\forCats\\" + testName);
            //ofstream out("D:\\6sem\\algorithm\\forCats\\" + output);

            in >> vertNum >> edgeNum >> fromVert >> endVert;
            startVert = fromVert;
            list.resize(vertNum);

            for (auto i = 0; i < vertNum; i++) {
                coordStruct mem{};
                in >> mem.x >> mem.y;
                coord.push_back(mem);
            }

            for (auto i = 0; i < edgeNum; i++) {
                int edgeStart, edgeEnd;
                float edgeWeight;
                in >> edgeStart >> edgeEnd >> edgeWeight;

                part mem{};
                mem.vert = edgeEnd;
                mem.weight = edgeWeight;

                list[edgeStart].push_back(mem);
            }

            start_time = high_resolution_clock::now();; // начальное время
            initialize();

            computeShortestPath();

            while (startVert != endVert) {

                if (dist[startVert].g == MAX) {//Наверное fromVert
                    break;
                }

                startVert = int(minSuccVert(startVert, 'v'));

                char choice;
                in >> choice;
                if ((choice == 'Y') && in.peek() != EOF) {
                    bool correctTest = true;
                    int changedVertNumber;
                    int changedEdgesNumber;
                    in >> changedVertNumber >> changedEdgesNumber;
                    vector<int> verts(0);


                    for (auto i = 0; i < changedVertNumber; i++) {
                        int vert, x, y;
                        in >> vert;// >> x >> y;
                        if (vert == startVert) {
                            /*coord[vert].x = x;
                            coord[vert].y = y;*/
                            verts.push_back(vert);
                            counter++;
                        }
                    }

                    for (auto i = 0; i < changedEdgesNumber; i++) {//correct
                        int start, end;
                        float changedWeight;
                        in >> start >> end >> changedWeight;

                        if (inVerts(verts, start, end)) {
                            //high_resolution_clock::time_point startMinusTime = high_resolution_clock::now();
                            for (auto j = 0; j < list[start].size(); j++) {
                                if (list[start][j].vert == end) {
                                    list[start][j].weight = changedWeight;
                                }
                            }
                            //high_resolution_clock::time_point endMinusTime = high_resolution_clock::now();
                            //duration_cast<duration<double>>(endMinusTime - startMinusTime);
                            updateVertex(start);
                        }
                    }

                    priority_queue<forQueue, vector<forQueue>, Compare> copyQ;
                    while (!Q.empty()) {
                        int qTopVert = Q.top().v;
                        Q.pop();
                        copyQ.push(calcKey(qTopVert));
                    }
                    Q = copyQ;
                    computeShortestPath();
                    //double res = Dijkstra();
                }

            }

            //unsigned int end_time = clock();
            end_time = high_resolution_clock::now();

            if (dist[fromVert].rhs == MAX) {
                dist[fromVert].rhs = -1;
            }

            duration<double> time = fullDijTime;
            duration<double> timeDNew = duration_cast<duration<double>>(end_time - start_time);
            if (abs(dist[fromVert].rhs - Dijkstra()) < eps) {
                cout << testNum << " " << " " << dist[fromVert].rhs << " is passed! D time: " << timeDNew.count() << "\n";
                     //out << dist[fromVert].rhs;
            } else {
                cout << " " << testNum << " !!!!!!!!!!!!Wrong answer!!!!!!!!!!!!!!! " << dist[fromVert].rhs << " "
                     << Dijkstra() << "\n";
            }
            times.push_back(timeDNew);
            //timeD += end_time - start_time;
            //out.close();
        }
        double t = 0;
        for (auto i: times) {
            t += i.count();
        }
        times.clear();
        cout << testNum << "Average time:" << t/10 << "\n";
    //}
    //cout << "Average time is: " << timeD/(testNum + 1) << "\n";

}
