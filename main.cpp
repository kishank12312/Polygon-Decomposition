/* Prithvi Hegde, 2020A7PS1718H,
Shashank Pandey, 2020A7PS0238H,
Kishan Koundinya, 2020A7PS0213H
*/

#include <bits/stdc++.h>
#include <limits.h>
using namespace std;

typedef double coordinate;

class halfEdge;
class vertex;
class face;
class DCEL;
vector<vertex> oginput;
/**
 * @brief Class definition for a Half Edge in the DCEL Data Structure 
*/
class halfEdge{
    public:
    vertex* origin;
    halfEdge* twin;
    face* incidentFace;
    halfEdge* nextHalfEdge;
    halfEdge* prevHalfEdge;
    float angle;

    /**
     * @brief Constructor of the halfEdge class
     * 
     * @param v1 One endpoint of the edge
     * @param v2 The other endpoint of the edge
     * 
     * @return void
    */
    halfEdge(vertex* v1, vertex* v2){
        origin = v1; 
        twin = NULL;
        incidentFace = NULL;
        nextHalfEdge = NULL;
        prevHalfEdge = NULL;
    }

};

/**
 * @brief Class definiton for a Vertex in the DCEL Data Structure
*/
class vertex{
    public:
    coordinate x, y;
    halfEdge* incidentEdge;
    int score;
    
    /**
     * @brief Constructor of the vertex class
     * 
     * @param a X coordinate
     * @param b Y coordinate
     * 
     * @return void
    */
    vertex(coordinate a, coordinate b){
        incidentEdge = NULL;
        x = a; y = b;
    }

    /**
     * @brief Default constructor of the vertex class
     * 
     * @return void
    */
    vertex(){}
    bool operator== (vertex const& v)const noexcept{
        return (v.x == x and v.y==y);
    }
    bool operator!= (vertex const& v)const noexcept{
        return !(v.x == x and v.y==y);
    }
    bool operator> (vertex const& v1)const{
        return (score > v1.score);
    }
 
    bool operator< (vertex const& v1) const{
        return (score < v1.score);
    }
 
    bool operator>= (vertex const& v1)const{
        return (score >= v1.score);
    }
 
    bool operator<= (vertex const& v1)const{
        return (score <= v1.score);
    }
};

/*
 * @brief Class definiton for a Face in the DCEL Data Structure
*/
class face{
    public:
    halfEdge* edge;
    bool external;
    
    /**
     * @brief Constructor of the face class
     * 
     * @return void
    */
    face(){
        edge = NULL;
        external = NULL;
    }

    /**
     * @brief Function that returns the vertices of the polygon represented by the face
     * 
     * 
     * @return vector<vertex*>
    */
    vector<vertex*> vertexList(){
        halfEdge* h = edge;
        vector<vertex*> ans;
        vertex* horig = (h->origin);
        ans.push_back(horig);

        while(h->nextHalfEdge != edge){
            h = h->nextHalfEdge;
            ans.push_back(h->origin);
        }
        return ans;
    }
};

/**
 * @brief Next Function: Takes a polygon and a vertex of the polygon, and returns the next vertex in the polygon
 * 
 * @param P It is a vector of vertex objects, denotes the vertices of the polygon in clockwise order
 * @param v It is a vertex of the polygon P, whose next vertex is required to be returned.
 * 
 * @return vertex
*/
vertex next(vector<vertex> P, vertex v){
    auto it = find(P.begin(), P.end(), v);
    it++;
    if(it == P.end()){
        return *P.begin();
    }
    else{
        return *it;
    }
}

/**
 * @brief Previous Function: Takes a polygon and a vertex of the polygon, and returns the previous vertex in the polygon
 * 
 * @param P It is a vector of vertex objects, denotes the vertices of the polygon in clockwise order
 * @param v It is a vertex of the polygon P, whose previous vertex is required to be returned.
 * 
 * @return vertex
*/
vertex prev(vector<vertex> P, vertex v){
    vertex fraud2;
    for(auto it = P.begin(); it != P.end(); ++it){
        if(*it == v){
            if (it == P.begin()) return *P.rbegin();
            else return *--it;
        }
    }
    return fraud2;
}

/**
 * @brief Class definition for the DCEL data structure used. 
*/
class DCEL{
    public:
    vector<vertex*> vertices;
    vector<halfEdge*> halfEdges;
    vector<face*> faces;
    
    /**
     * @brief Defualt constructor of the DCEL class
     * 
     * @return void
    */
    DCEL(){}

    /**
     * @brief Function that builds a DCEL of a polygon with no diagnals
     * 
     * @param vl List of vertices of the polygon in order
     * 
     * @return void
    */
    void buildDCEL(vector<vertex> vl){
        //Creation of vertexList
        for(auto i: vl){
            vertex* v = new vertex(i.x, i.y);
            v->score = i.score;
            vertices.push_back(v);
        }

        //Creation of halfEdgelist + assignment of twin
        for (int i = 0; i < vl.size(); i++)
        {
            vertex* nvert;
            if(i == vl.size()-1)
                nvert = vertices[0];
            else
                nvert = vertices[i+1];
            halfEdge* h1 = new halfEdge(vertices[i], nvert);
            halfEdge* h2 = new halfEdge(nvert, vertices[i]);
            h1->twin = h2;
            h2->twin = h1;

            nvert->incidentEdge = h2;

            halfEdges.push_back(h2);
            halfEdges.push_back(h1);
        }


        //id of next and prev halfEdges
        for (int i = 0; i < vertices.size(); i++)
        {
            vertex* prevVertex;
            if(i==0)
                prevVertex = vertices[vertices.size() - 1];
            else
                prevVertex = vertices[i-1];
            vertices[i]->incidentEdge->nextHalfEdge = prevVertex->incidentEdge;
            
            vertex* nextvertex;
            if(i == vertices.size()-1)
                nextvertex = vertices[0];
            else
                nextvertex = vertices[i+1];
            vertices[i]->incidentEdge->prevHalfEdge = nextvertex->incidentEdge;

            vertices[i]->incidentEdge->twin->nextHalfEdge = nextvertex->incidentEdge->twin;
            vertices[i]->incidentEdge->twin->prevHalfEdge = prevVertex->incidentEdge->twin;



        }
        

        //4. Assigning faces and check if external
        face* f1 = new face();
        face* f2 = new face();

        f1->edge = vertices[0]->incidentEdge;
        f2->edge = vertices[0]->incidentEdge->twin;

        f1->external = false;
        f2->external = true;
        faces.push_back(f1);
        faces.push_back(f2);
        return;
    }

    //5. Add edge between two vertices
    /**
     * @brief Function that adds an edge in the DCEL between 2 vertices
     * 
     * @param v1 One of the endpoints of the edge
     * @param v2 The other endpoint of the edge
     * 
     * @return void
    */
    void addEdge(vertex* v1, vertex* v2){
        //check if edge already present
        for(auto e: halfEdges){
            if (*(e->origin) == *v1 and *(e->nextHalfEdge->origin) == *v2)
                return;
        }

        int targetFace = -1;
        for (int i = 0; i < faces.size(); i++)
        {
            face* f = faces[i];
            if(f->external) continue;

            vector<vertex*> vList = f->vertexList();
            int count = 0;
            for(auto v: vList){
                if(*v == *v1){
                    v1 = v;
                    count++;
                }
                if(*v == *v2){
                    v2 = v;
                    count++;
                }
            }
            if(count==2){
                targetFace = i;
                break;
            }
        }
        if(targetFace == -1){
            perror("Target face not found");
            exit(1);
        }

        halfEdge* e1 = faces[targetFace]->edge;
        halfEdge* e2 = faces[targetFace]->edge;
        halfEdge* temp = faces[targetFace]->edge;
        do{
            if(temp->twin->origin == v1)
                e1 = temp;
            else if(temp->origin == v2)
                e2 = temp;
            temp = temp->nextHalfEdge;
        } while(temp != faces[targetFace]->edge);
        halfEdge *e3 = e1->nextHalfEdge;
        halfEdge *e4 = e2->prevHalfEdge;

        //create new edges and update next nad prev
        halfEdge* h1 = new halfEdge(v1, v2); 
        halfEdge* h2 = new halfEdge(v2, v1); 

        h1->twin = h2;
        h2->twin = h1;

        h1->nextHalfEdge = e2;
        e2->prevHalfEdge = h1;
        h2->nextHalfEdge = e3;
        e3->prevHalfEdge = h2;

        h1->prevHalfEdge = e1;
        e1->nextHalfEdge = h1;
        h2->prevHalfEdge = e4;
        e4->nextHalfEdge = h2;

        face* fnew = new face();
        face* fnew2 = new face();

        fnew->edge = h1;
        fnew2->edge = h2;

        fnew->external = false;
        fnew2->external = false;

        faces[targetFace] = fnew;
        faces.push_back(fnew2);
        halfEdges.push_back(h1);
        halfEdges.push_back(h2);
        return;
    }

    /**
     * @brief Function that removes an edge in the DCEL between 2 vertices
     * 
     * @param v1 One of the endpoints of the edge
     * @param v2 The other endpoint of the edge
     * @param originalpoly The original polygon
     * 
     * @return void
    */
    void removeEdge(vertex* v1, vertex* v2,vector<vertex> originalpoly){
        DCEL *newdcel = new DCEL();
        newdcel->buildDCEL(originalpoly);
        for(auto e: halfEdges){
            if((*(e->origin) == *(v1)) and (*(e->twin->origin) == *(v2))){
                continue;
            }
            else if((*(e->origin) == *(v2)) and (*(e->twin->origin) == *(v1))){
                continue;
            }
            else{
                newdcel->addEdge(e->origin,e->twin->origin);
            }
        }
        this->faces = newdcel->faces;
        this->halfEdges = newdcel->halfEdges;
        this->vertices = newdcel->vertices;
        //delete newdcel;
    }

    /**
     * @brief Function that finds the faces that share a diagnal with given endpoints
     * 
     * @param v1 One of the endpoints of the diagnal
     * @param v2 The other endpoint of the diagnal
     * 
     * @return vector<face>
    */
    vector<face> findfaces(vertex* v1, vertex* v2){
    //check if edge already present
        vector<int> targetFace;
        for (int i = 0; i < faces.size(); i++)
        {
            face* f = faces[i];
            if(f->external) continue;

            vector<vertex*> vList = f->vertexList();
            int count = 0;
            for(auto v: vList){
                if(*v == *v1){
                    v1 = v;
                    count++;
                }
                if(*v == *v2){
                    v2 = v;
                    count++;
                }
            }
            if(count==2){
                targetFace.push_back(i);
            }
        }
        if(targetFace.size()!=2){
            perror("Unable to find faces");
            exit(1);
        }
        else {
            vector<face> ans;
            for(auto i:targetFace){
                ans.push_back(*(faces[i]));
            }
            return ans;
        }

    }
    
    /**
     * @brief Function that finds the face to which 2 vertices belong
     * 
     * @param v1 One of the vertices
     * @param v2 The other vertex
     * 
     * @return face
    */
    face findface(vertex* v1, vertex* v2){
        int targetFace = -1;
        for (int i = 0; i < faces.size(); i++)
        {
            face* f = faces[i];
            if(f->external) continue;

            vector<vertex*> vList = f->vertexList();
            int count = 0;
            for(auto v: vList){
                if(*v == *v1){
                    v1 = v;
                    count++;
                }
                if(*v == *v2){
                    v2 = v;
                    count++;
                }
            }
            if(count==2){
                targetFace = i;
                break;
            }
        }
        if(targetFace == -1){
            perror("Target face not found");
            exit(1);
        }
        else return *(faces[targetFace]);
    }

    /**
     * @brief Function that outputs the DCEL's polygons into a file
     * 
     * @param filename Name of output file
     * 
     * @return void
    */
    void printDCEL(string filename){
        FILE* fp = fopen(filename.c_str(), "w");
        fprintf(fp, "%d\n", faces.size()-1);

        for(auto f: faces){
            if(f->external)continue;
            fprintf(fp, "\n");

            vector<vertex*> v2 = f->vertexList();
            fprintf(fp, "%d\n", v2.size());
            for(auto v: v2){
                fprintf(fp, "%f %f\n", v->x, v->y);
            }
        }
        fclose(fp);

        return;
    }
};

/**
 * @brief Function that takes three consecutive points which make an angle, returns whether that angle is reflex or not
 * 
 * @param a The first vertex of the angle
 * @param b The second vertex of the angle
 * @param c The third vertex of the angle
 * @return bool
*/
bool isNotReflex(vertex a, vertex b, vertex c){
    float res = (b.x - a.x) * (c.y - b.y) - (c.x - b.x) * (b.y - a.y);
    if(abs(res) < 1e-7){
        res = 0;
    }
    if(res <= 0)
        return true;
    return false;
}

/**
 * @brief Function that takes a ploygon and returns the smallest rectangle that contains the polygon
 * 
 * @param L The list of vertices of the polygon
 * @return vector<vertex>
*/
vector<vertex> rect(vector<vertex> L){
    coordinate Xmax=(coordinate)INT32_MIN, Xmin=(coordinate)INT32_MAX, Ymax=(coordinate)INT32_MIN, Ymin=(coordinate)INT32_MAX;
    for(auto i: L){
        if(i.x > Xmax){
            Xmax = i.x;
        }
        if(i.x < Xmin){
            Xmin = i.x;
        }
        if(i.y > Ymax){
            Ymax = i.y;
        }

        if(i.y < Ymin){
            Ymin = i.y;
        }
    }
    vertex v1(Xmin, Ymax), v2(Xmax, Ymax), v3(Xmax, Ymin), v4(Xmin, Ymin);
    vector<vertex> ans = {v1,v2,v3,v4};
    return ans;
}

/**
 * @brief Function that takes a rectangle, and a vertex and returns true if the vertex is in the rectangle
 * 
 * @param R The list of vertices of the rectangle
 * @param v The vertex
 * @return bool
*/
bool inRect(vector<vertex> R, vertex v){
    if((v.x <= R[1].x and v.x >= R[3].x) and (v.y <= R[1].y and v.y >= R[3].y)){
        return true;
    }
    return false;
}

/**
 * @brief Function that takes a polygon, and a vertex and returns true if the vertex is in the polygon
 * 
 * @param point The vertex
 * @param polygon The list of vertices of the polygon
 * @return bool
*/
bool point_in_polygon(vertex point, vector<vertex> polygon) {
   int num_vertices = polygon.size();
   double x = point.x, y = point.y;
   bool inside = false;
   vertex p1 = polygon[0], p2;
   for (int i = 1; i <= num_vertices; i++) {
       p2 = polygon[i % num_vertices];
       if (y >= min(p1.y, p2.y)) {
           if (y <= max(p1.y, p2.y)) {
               if (x <= max(p1.x, p2.x)) {
                   double x_intersection = (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
                   if (p1.x == p2.x || x <= x_intersection) {
                       inside = !inside;
                   }
               }
           }
       }
       p1 = p2;
   }
   return inside;
}

/**
 * @brief Function that takes 3 points, and returns on which side the 3rd point lies of the line made by first 2 points
 * 
 * @param p1 One of the endpoints of the line
 * @param p2 The other endpoint of the line
 * @param v The vertex whose location is to be tested
 * @return bool
*/
bool whichSide(vertex p1, vertex p2, vertex v){
    float res = ((v.x - p1.x)*(p2.y - p1.y)-(v.y - p1.y)*(p2.x - p1.x));
    if(abs(res) < 1e-7){
        res =0;
    }
    return 0 > res ;
}

/**
 * @brief Function that returns the angle made by the line joining two coordinate pairs(Points)
 * 
 * @param center One of the endpoints of the line
 * @param p The other endpoint of the line
 * 
 * @return float
*/
float getAngle(vector<coordinate> center, vector<coordinate> p){
    float x = p[0] - center[0];
    float y = p[1] - center[1];
    float angle = atan2(y,x);
    angle = angle * 180 / 3.141592;
    if(angle <= 0){
        return angle + 360;
    }
    else
        return angle;
}

/**
 * @brief Function that returns the distance between two coordinate pairs(Points)
 * 
 * @param p1 One of the points
 * @param p2 The other point
 * 
 * @return float
*/
float getDistance(vector<coordinate> p1, vector<coordinate> p2){
    float x = p1[0] - p2[0];
    float y = p1[1] - p2[1];
    return x*x + y*y;
}

/**
 * @brief Struct that is used as the comparator for the sortPoints fucntions
 * 
*/
struct comparePoints{
    pair<float,float> center;
    comparePoints(pair<float,float> center){
        this->center = center;  
    }
    inline bool operator() (const vector<coordinate>& p1, const vector<coordinate>& p2)
    {
        float angle1 = getAngle({0,0},p1);
        float angle2 = getAngle({0,0},p2);
        if(angle1<angle2){
            return true;
        }
        float d1 = getDistance({0,0},p1);
        float d2 = getDistance({0,0},p2);
        if((angle1==angle2) and (d1<d2)){
            return true;
        }
        return false;
    }
};

/**
 * @brief Function that sorts the points of a given polygon in clockwise order
 * 
 * @param points Points of the polygon
 * 
 * @return vector<vector<coordinate>>
*/
vector<vector<coordinate>> sortPoints(vector<vector<coordinate>> points){
    float centroidx=0, centroidy = 0;
    for(auto i:points){
        centroidx += i[0];
        centroidy += i[1];
    }
    centroidx = centroidx/points.size();
    centroidy = centroidy/points.size();
    for (int i = 0; i < points.size(); i++)
    {
        points[i][0] -= centroidx;
        points[i][1] -= centroidy;
    }
    sort(points.begin(), points.end(),comparePoints({centroidx,centroidy}));
    for (int i = 0; i < points.size(); i++)
    {
        points[i][0] += centroidx;
        points[i][1] += centroidy;
    }
    return points;
    
}

/**
 * @brief Function that implements algorithm 1, decomposition of a polygon into convex polygons, stored in a DCEL
 * 
 * @param P Points of the polygon in clockwise order
 * 
 * @return DCEL*
*/
DCEL* algo1(vector<vertex> P){
    DCEL *ans = new DCEL();
    ans->buildDCEL(P);

    // 2
    vector<vector<vertex>> allL;
    vector<vertex> L;
    L.push_back(*P.begin());
    allL.push_back(L);
    int m=1;
    int n=P.size();


    // 3
    while(n>3){
        // 3.1
        vertex v1=*(allL[m-1].rbegin());
        vertex v2 = next(P, v1);

        // 3.2
        vector<vertex> vArray = {vertex(),v1,v2};
        allL.push_back({v1, v2});
        int i=2;
        vArray.push_back(next(P, vArray[i]));

        // 3.3
        while((allL[m].size()<n) 
            && (isNotReflex(vArray[i-1], vArray[i], vArray[i+1]) )
            && (isNotReflex(vArray[i], vArray[i+1], vArray[1])) 
            && (isNotReflex(vArray[i+1], vArray[1], vArray[2]))
            )
            {

            // 3.3.1
            allL[m].push_back(vArray[i+1]);
            i++;
            vArray.push_back(next(P, vArray[i]));
        }


        // 3.4
        if(allL[m].size() != P.size()){

            // 3.4.1
            vector<vertex> LPVS;
            vector<vertex> allL2 = allL[m];
            sort(allL2.begin(), allL2.end()); 
            set_difference(begin(P), end(P), begin(allL2), end(allL2), inserter(LPVS, LPVS.end()));
            vector<vertex> realLPVS;
            for(auto i:LPVS){
                vertex p = prev(oginput,i);
                vertex n = next(oginput,i);
                if(isNotReflex(p,i,n))
                    continue;
                else
                    realLPVS.push_back(i);
            }
            LPVS = realLPVS;

            // 3.4.2
            while(LPVS.size() > 0){

                vector<vertex> R = rect(allL[m]);

                bool backward = false;

                while(!backward and (LPVS.size()>0)){
                    vertex v;
                    do{
                        v = *LPVS.begin();

                        if(!inRect(R, v)){
                            LPVS.erase(find(LPVS.begin(), LPVS.end(), v));
                        }
                    }while(!inRect(R, v) and (LPVS.size()!=0));

                    if(LPVS.size()>0){

                        if(point_in_polygon(v,allL[m])){
                            bool side = whichSide(vArray[1], v, *allL[m].rbegin());
                            vector<vertex> VTR;
                            for(auto i: allL[m]){
                                if(side == whichSide(vArray[1], v, i)){
                                    VTR.push_back(i);
                                }
                            }
                            vector<vertex> newLm;
                            vector<vertex> allL2 = allL[m];
                            vector<vertex> VTR2 = VTR;
                            sort(allL2.begin(), allL2.end()); 
                            sort(VTR2.begin(), VTR2.end()); 
                            set_difference(begin(allL2), end(allL2), begin(VTR2), end(VTR2), inserter(newLm, newLm.end()));
                            allL[m] = newLm;
                            backward = true;
                        }
                    LPVS.erase(find(LPVS.begin(), LPVS.end(), v));
                    }
                }
            }
        }
        
        // 3.5
        vertex lastlLm = *allL[m].rbegin();
        if(lastlLm != vArray[2]){
            
            // 3.5.1 ans
            for (int i = 0; i < allL[m].size(); i++){
                int n = i+1;
                if(n==allL[m].size())
                    n = 0;
                ans->addEdge(&allL[m][i],&allL[m][n]);
            }
            
            // 3.5.2
            vertex first = *begin(allL[m]);
            vertex last = *(allL[m].rbegin());
            vector<vertex> newPP;
            vector<vertex> allL2 = allL[m];
            sort(allL2.begin(), allL2.end());
            set_difference(begin(P), end(P), begin(allL2), end(allL2), inserter(newPP, newPP.end()));
            P = newPP;
            P.push_back(first);
            P.push_back(last);
            set<vertex> tempset(P.begin(),P.end());
            vector<vertex> awe(tempset.begin(),tempset.end());
            P = awe;
            sort(P.begin(), P.end());
            n = n - (allL[m].size()) + 2;

        }

        // 3.6
        m++;
    }
    return ans;
}

/**
 * @brief Function that takes a DCEL and returns all the internal faces in that DCEL represented as vector of vertices
 * 
 * @param dcel The input DCEL
 * 
 * @return vector<vector<vertex>>
*/
vector<vector<vertex>> polygonList(DCEL dcel){
    vector<vector<vertex>> pList;

    for(auto f: dcel.faces){
        if(f->external) continue;
        vector<vertex*> vListPtr = f->vertexList();
        vector<vertex> vList; 
        for(auto v:vListPtr){
            vList.push_back(*v);
        }
        reverse(vList.begin(),vList.end());
        pList.push_back(vList);

    }
    return pList;
}

/**
 * @brief Function that checks if the edge made by two vertices is a diagnal in the original polygon
 * 
 * @param v1 One of the endpoints of the edge
 * @param v2 The other endpoint of the edge
 * 
 * @return bool
*/
bool isDiagonal(vertex v1, vertex v2){
    if(v2 == next(oginput, v1)){
        return false;
    }
    else if(v2 == prev(oginput, v1)){
        return false;
    }
    else if(v1 == next(oginput, v2)){
        return false;
    }
    else if(v1 == prev(oginput, v2)){
        return false;
    }
    else{
        return true;
    }
}

/**
 * @brief Function that checks if diagnal made by two vertices is in the given polynomial
 * 
 * @param polygon The polygon
 * @param v1 One of the endpoints of the diagnal
 * @param v2 The other endpoint of the diagnal
 * 
 * @return bool
*/
bool isInPolygon(vector<vertex> polygon, vertex v1, vertex v2){
    if(find(polygon.begin(), polygon.end(), v1) != polygon.end()){
        if(find(polygon.begin(), polygon.end(), v2) != polygon.end()){
            return true;
        }
    }
    return false;
}

/**
 * @brief Function that generates the list LLE as described in the paper
 * 
 * @param polygon List of polygons in the DCEL
 * @param dcel The DCEL of the decomposed polygon
 * 
 * @return map<pair<vertex, vertex>, vector<int>>
*/
map<pair<vertex, vertex>, vector<int>> LLE(vector<vector<vertex>> polygon, DCEL dcel){
    map<pair<vertex, vertex>, vector<int>> lle;
    for(auto e: dcel.halfEdges){
        vertex v1 = *(e->origin);
        vertex v2 = *(e->twin->origin);

        if(lle.count({v2, v1})){
            continue;
        }

        //check if diag
        if(!isDiagonal(v1, v2)){
            continue;
        }
        else{
            lle[{v1, v2}] = {};

            for(int i = 0; i < polygon.size(); i++){
                vector<vertex> p = polygon[i];
                if(isInPolygon(p, v1, v2)){
                    lle[{v1, v2}].push_back(i);
                }
            }
        }
    }
    return lle;
} 

/**
 * @brief Function that generates the list LP as described in the paper
 * 
 * @param lle The list LLE as described in the paper
 * @return map<vertex, vector<pair<int, vertex>>>
*/
map<vertex, vector<pair<int, vertex>>> LP(map<pair<vertex, vertex>,vector<int>> lle){
    map<vertex, vector<pair<int, vertex>>> lp;
    for(auto i:oginput){
        lp[i] = {};
    }
    for(auto l: lle){
        vertex v1, v2;
        v1 = l.first.first;
        v2 = l.first.second;

        for(auto j: l.second){
            lp[v1].push_back({j, v2});
            lp[v2].push_back({j, v1});
        }
    }
    return lp;
}

/**
 * @brief Function that determines whether a given vertex is not a notch
 * 
 * @param v The vertex
 * @return bool
*/
bool isConvex(vertex v){
    vertex p = prev(oginput,v);
    vertex n = next(oginput,v);
    if(isNotReflex(p,v,n))
        return true;
    else
        return false;
}

/**
 * @brief Function joins 2 polygons separated by a diagnal by removing the diagnal
 * 
 * @param dcel DCEL of the polygons
 * @param vs One of the endpoints of the diagnal
 * @param vt The other endpoint of the diagnal
 * @return vector<vertex>
*/
vector<vertex> joinPoly(DCEL dcel, vertex vs, vertex vt){
    face f = dcel.findface(&vs,&vt);

    vector<vertex*> vListPtr = f.vertexList();
    vector<vertex> vList; 
    for(auto v:vListPtr){
        vList.push_back(*v);
    }
    reverse(vList.begin(),vList.end());
    return vList;
}

//Merging Process
/**
 * @brief Function implements the Merging Process Algorithm in the paper
 * 
 * @param dcel DCEL of the polygons
 * @param lle The list lle as described in the paper
 * @param lp The list lp as described in the paper
 * 
 * @return DCEL
*/
DCEL merging(DCEL dcel,vector<vector<vertex>> poly, map<pair<vertex, vertex>, vector<int>> lle, map<vertex, vector<pair<int, vertex>>> lp){
    vector<bool> LDP;
    vector<int> LUP;
    int NP, m = lle.size();

    // 1.
    NP = m+1;

    // 2.
    for (int i = 0; i <= NP; i++){
        LDP.push_back(true);
        LUP.push_back(i);
        
    }

    // 3.
    for (auto kvp: lle)
    {
        // 3.1

        vertex vs, vt;
        vs = kvp.first.first;
        vt = kvp.first.second;

        vector<int> joiningPolygons = kvp.second;
        //3.2
        if((lp[vs].size() > 2 and lp[vt].size()>2)
        or (lp[vs].size()>2 and isConvex(vt)) 
        or (lp[vt].size()>2 and isConvex(vs)) 
        or (isConvex(vs) and isConvex(vt))){

            //3.2.1
            vertex j1, j2, j3, i1, i2, i3;
            j2 = vt;
            i2 = vs;
            vector<face> faces = dcel.findfaces(&j2, &i2);
            halfEdge* temp1 = faces[0].edge;
            do{
                if(*(temp1->origin) == vs){
                    if(*(temp1->twin->origin) != vt){
                        i1 = *(temp1->twin->origin);
                    }
                }
                if(*(temp1->twin->origin) == vt){
                    if(*(temp1->origin) != vs){
                        j3 = *(temp1->origin);
                    }
                }
                if(*(temp1->twin->origin) == vs){
                    if(*(temp1->origin) != vt){
                        i3 = *(temp1->origin);
                    }
                }
                if(*(temp1->origin) == vt){
                    if(*(temp1->twin->origin) != vs){
                        j1 = *(temp1->twin->origin);
                    }
                }
                temp1 = temp1->nextHalfEdge;

            } while(temp1 != faces[0].edge);
            temp1 = faces[1].edge;
            do{
                if(*(temp1->origin) == vs){
                    if(*(temp1->twin->origin) != vt){
                        i1 = *(temp1->twin->origin);
                    }
                }
                if(*(temp1->twin->origin) == vt){
                    if(*(temp1->origin) != vs){
                        j3 = *(temp1->origin);
                    }
                }
                if(*(temp1->twin->origin) == vs){
                    if(*(temp1->origin) != vt){
                        i3 = *(temp1->origin);
                    }
                }
                if(*(temp1->origin) == vt){
                    if(*(temp1->twin->origin) != vs){
                        j1 = *(temp1->twin->origin);
                    }
                }
                temp1 = temp1->nextHalfEdge;

            } while(temp1 != faces[1].edge);


            //3.2.2
            //3.2.3
            //3.2.4
            if(isNotReflex(i1, i2, i3) and (isNotReflex(j1,j2,j3))){

                NP++;

                // Write
                dcel.removeEdge(&vs,&vt, oginput);
                vector<vertex> newp = joinPoly(dcel,vs,vt);
                poly.push_back(newp);

                // Updating LDP and LUP

                LUP[LUP[joiningPolygons[0]]] = NP;
                LUP[LUP[joiningPolygons[1]]] = NP;

                for (int h = 0; h < poly.size(); h++)
                {
                    if((LUP[h] == joiningPolygons[0]) or (LUP[h] == joiningPolygons[1])){
                        LUP[h] = NP;
                    }
                }
            }
        }
    }
    dcel.printDCEL("Merged Polygon.txt");
    exit(1);
    return dcel;
}

/**
 * @brief Main function
 * 
 * @param argc Number of command line arguments provided
 * @param argv List of the command line arguments provided
 * 
 * @return int
*/
int main(int argc, char* argv[]){
    FILE *fp = freopen(argv[1], "r", stdin);
    vector<vector<coordinate>> dataPoints;

    string buffer;
    int edgesCount;
    cin >> edgesCount;
    // edgesCount--;
    while(edgesCount){
        vector<coordinate> coords(2);
        cin >> coords[0] >> coords[1];
        edgesCount--;
        dataPoints.push_back(coords);
    }
    // dataPoints = sortPoints(dataPoints);
    // reverse(dataPoints.begin(),dataPoints.end());

    vector<vertex> input;
    int score = 1;
    for(auto i: dataPoints){
        vertex v(i[0],i[1]);
        v.score = score;
        score++;
        input.push_back(v);
    }
    oginput = input;
    DCEL *output = algo1(input);
    
    output->printDCEL("Unmerged Polygon.txt");

    //MERGING ALGO
    vector<vector<vertex>> poly = polygonList(*output);
    map<pair<vertex, vertex>, vector<int>> lle = LLE(poly, *output);
    map<vertex, vector<pair<int, vertex>>> lp = LP(lle);
    DCEL lup = merging(*output, poly, lle, lp);

    lup.printDCEL("Merged Polygon.txt");
    
    return 0;
}

