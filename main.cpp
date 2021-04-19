#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <set>

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

int vertNum, edgeNum, startVert, endVert;

priority_queue<forQueue, vector<forQueue>, Compare> Q;
vector <vector <part>> list(0);
vector <distStruct> dist (0);
vector <coordStruct> coord(0);
int fromVert = 0;
const auto MAX = 100000;
const auto eps = 0.001;



double Dijkstra() {
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

int operator< (forQueue a, forQueue b) { //tested
    if (a.key1 < b.key1) {
        return true;
    }
    if (a.key1 > b.key1) {
        return false;
    }
    if (a.key2 < b.key2) {
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
        pathToS = getEdgeValue(prevVert, s) + dist[s].rhs;
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
        distStruct mem{};
        mem.g = MAX;
        mem.rhs = MAX;
        dist.push_back(mem);
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
        (dist[fromVert].rhs != dist[fromVert].g)) && (!Q.empty()))   {
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

    //for (int i = 0; i < 1000; i++) {
        clearData();
        int testNum = 332;
        string testName = "input" + to_string(testNum) + ".txt";
        ifstream in("D:\\6sem\\algorithm\\tests4\\" + testName);
        ofstream out("output.txt");

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

        /*for (auto i = 0; i < list.size(); i++){
            for (auto j: list[i]) {
                cout << i << " " << j.vert << " "<< j.weight<< "  ";
            }
            cout << coord[i].x << " " << coord[i].y << "     ";
        }*/

        initialize();

        computeShortestPath();

        /*forQueue mem{};
        mem.v = 1;
        mem.key1 = 4;
        mem.key2 = 3;
        Q.push(mem);

        forQueue mem2{};
        mem2.v = 2;
        mem2.key1 = 5;
        mem2.key2 = 8;
        Q.push(mem2);

        mem2.v = 3;
        mem2.key1 = 4;
        mem2.key2 = 4;
        Q.push(mem2);

        cout << "   ";

        remove(2);
        while (!Q.empty()) {
            cout << Q.top().v << " ";
            Q.pop();
        }*/

        while (startVert != endVert) {

            //cout << startVert << "\n" << " " << dist[0].rhs << " " << dist[0].g << "\n" << dist[1].rhs << " " << dist[1].g << "\n"<< dist[2].rhs << " " << dist[2].g << "\n"<< dist[3].rhs << " " << dist[3].g << "\n\n";
            if (dist[startVert].g == MAX) {
                break;
            }

            startVert = int(minSuccVert(startVert, 'v'));

            char choice;
            in >> choice;
            if ((choice == 'Y') && in.peek() != EOF) {
                int changedVertNumber;
                int changedEdgesNumber;
                in >> changedVertNumber >> changedEdgesNumber;
                for (auto i = 0; i < changedVertNumber; i++) {
                    int vert, x, y;
                    in >> vert >> x >> y;
                    coord[vert].x = x;
                    coord[vert].y = y;
                }

                /*COUT FOR COORDINATES
                 for (auto i: coord) {
                    cout << "\n" << i.x << " " << i.y;
                }*/

                for (auto i = 0; i < changedEdgesNumber; i++) {//correct
                    int start, end;
                    float changedWeight;
                    in >> start >> end >> changedWeight;
                    for (auto j = 0; j < list[start].size(); j++) {
                        if (list[start][j].vert == end) {
                            list[start][j].weight = changedWeight;
                        }
                    }
                    updateVertex(start);
                }

                /*COUT FOR LIST
                 for (int i = 0; i < list.size(); i++) {
                    for (auto j: list[i]) {
                        cout << "\n" << i << " " << j.vert << " " << j.weight ;
                    }
                }
                cout << "\n\n";*/

             /*CHECK IF QUEUE UPDATES VALUE OF CALK KEY AND IT DOESNT
              priority_queue<forQueue, vector<forQueue>, Compare> copyQTest;
            copyQTest = Q;
            while (!copyQTest.empty()) {
                cout << "\n" << copyQTest.top().v << " " << copyQTest.top().key1 << " " << copyQTest.top().key2;
                copyQTest.pop();
            }
            copyQTest = Q;
            priority_queue<forQueue, vector<forQueue>, Compare> copyQTest2;
            while (!copyQTest.empty()) {
                int qTopVert = copyQTest.top().v;
                copyQTest.pop();
                copyQTest2.push(calcKey(qTopVert));
            }
            while (!copyQTest2.empty()) {
                cout << "\n" << copyQTest2.top().v << " " << copyQTest2.top().key1 << " " << copyQTest2.top().key2;
                copyQTest2.pop();
            }*/

                priority_queue<forQueue, vector<forQueue>, Compare> copyQ;
                while (!Q.empty()) {
                    int qTopVert = Q.top().v;
                    Q.pop();
                    copyQ.push(calcKey(qTopVert));
                }
                Q = copyQ;
                computeShortestPath();
            }
        }

        /*for (auto i: dist) {
            cout << "\n" << i.rhs << " " << i.g;
        }*/

        if (dist[fromVert].rhs == MAX) {
            dist[fromVert].rhs = -1;
        }

        if (abs(dist[fromVert].rhs - Dijkstra()) < eps) {
            //cout << "Passed!\n";
        } else {
            cout << "!!!!!!!!!!!!Wrong answer!!!!!!!!!!!!!!! " << dist[fromVert].rhs << " " << Dijkstra() << " " << testNum << "\n";
        }
        //cout << "   " << dist[fromVert].rhs << " " << Dijkstra() << "\n";


        /*vector<int> res = pred(5);
        for (auto i : res) {
            cout << i << " ";
        }

        vector<int> res = succ(5);
        for (auto i : res) {
            cout << i << " ";
        }*/
    //}

}
