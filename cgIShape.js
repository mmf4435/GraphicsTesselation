//
// fill in code that creates the triangles for a cube with dimensions 1x1x1
// on each side (and the origin in the center of the cube). with an equal
// number of subdivisions along each cube face as given by the parameter
//subdivisions
//
function makeCube (subdivisions)  {
    //left bottom front
    let p1 = [-0.5, -0.5, 0.5];
    //right bottom front
    let p2 = [0.5, -0.5, 0.5];
    //left top front
    let p3 = [-0.5, 0.5, 0.5];
    //right top front
    let p4 = [0.5, 0.5, 0.5];
    //left bottom back
    let p5 = [-0.5, -0.5, -0.5];
    //right bottom back
    let p6 = [0.5, -0.5, -0.5];
    //left top back
    let p7 = [-0.5, 0.5, -0.5];
    //right top back
    let p8 = [0.5, 0.5, -0.5];

    tessSquare(subdivisions, p1, p2, p3, p4);
    tessSquare(subdivisions, p6, p5, p8, p7);

    tessSquare(subdivisions, p2, p6, p4, p8);
    tessSquare(subdivisions, p5, p1, p7, p3);

    tessSquare(subdivisions, p3, p4, p7, p8);
    tessSquare(subdivisions, p5, p6, p1, p2);
}

//recursively subdivide and create triangles
function tessSquare(subdivisions, pt1, pt2, pt3, pt4){
    if(subdivisions == 1){
        addTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
        addTriangle(pt2[0], pt2[1], pt2[2], pt4[0], pt4[1], pt4[2], pt3[0], pt3[1], pt3[2]);
    }
    else{
        let q12 = calcMidpoint(pt1, pt2);
        let q13 = calcMidpoint(pt1, pt3);
        let q24 = calcMidpoint(pt2, pt4);
        let q34 = calcMidpoint(pt3, pt4);
        let qCenter = calcMidpoint(pt1, pt4);
        tessSquare(subdivisions - 1, pt1, q12, q13, qCenter);
        tessSquare(subdivisions - 1, q12, pt2, qCenter, q24);
        tessSquare(subdivisions - 1, q13, qCenter, pt3, q34);
        tessSquare(subdivisions - 1, qCenter, q24, q34, pt4);
    }
}

function calcMidpoint(pt1, pt2){
    return [(pt1[0] + pt2[0])/2, (pt1[1] + pt2[1])/2, (pt1[2] + pt2[2])/2];
}


//
// fill in code that creates the triangles for a cylinder with diameter 1
// and height of 1 (centered at the origin) with the number of subdivisions
// around the base and top of the cylinder (given by radialdivision) and
// the number of subdivisions along the surface of the cylinder given by
//heightdivision.
//
function makeCylinder (radialdivision,heightdivision){
    //radial setup
    let r = 0.5;
    let ptCenterTop = [0, 0.5, 0];
    let ptCenterBottom = [0, -0.5, 0];
    let divInc = 360 / radialdivision;
    let triSplit = divInc / 2;
    let alpha = 0;

    let pt2Top = [r*Math.cos(radians(alpha)), 0.5, r*Math.sin(radians(alpha))];
    let pt2Bottom = [r*Math.cos(radians(alpha)), -0.5, r*Math.sin(radians(alpha))];
    while(alpha <= 360){
        //sides setup
        let vertDiff = ptCenterTop[1] - ptCenterBottom[1];
        let vertInc = vertDiff / heightdivision;
        let vertLoc = ptCenterBottom[1];

        //handle radial
        alpha += triSplit;
        let pt1Top = [r*Math.cos(radians(alpha)), 0.5, r*Math.sin(radians(alpha))];
        let pt1Bottom = [r*Math.cos(radians(alpha)), -0.5, r*Math.sin(radians(alpha))];
        addTriangle(ptCenterTop[0],ptCenterTop[1],ptCenterTop[2], pt1Top[0], pt1Top[1], pt1Top[2], pt2Top[0], pt2Top[1], pt2Top[2]);
        addTriangle(ptCenterBottom[0],ptCenterBottom[1],ptCenterBottom[2], pt2Bottom[0], pt2Bottom[1], pt2Bottom[2], pt1Bottom[0], pt1Bottom[1], pt1Bottom[2]);

        //handle sides
        let pt1 = pt1Bottom;
        let pt2 = pt2Bottom;
        while(vertLoc <= ptCenterTop[1]){
            let pt3 = [pt1[0], vertLoc, pt1[2]];
            let pt4 = [pt2[0], vertLoc, pt2[2]];
            addTriangle(pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], pt3[0], pt3[1], pt3[2]);
            addTriangle(pt2[0], pt2[1], pt2[2], pt4[0], pt4[1], pt4[2], pt3[0], pt3[1], pt3[2]);
            pt1 = pt3;
            pt2 = pt4;
            vertLoc += vertInc;
        }
        pt2Top = pt1Top;
        pt2Bottom = pt1Bottom;
    }
}


//
// fill in code that creates the triangles for a cone with diameter 1
// and height of 1 (centered at the origin) with the number of
// subdivisions around the base of the cone (given by radialdivision)
// and the number of subdivisions along the surface of the cone
//given by heightdivision.
//
function makeCone (radialdivision, heightdivision) {
    let ptApex = [0, 0.5, 0];

    //radial setup
    let r = 0.5;
    let ptCenterRad = [0, -0.5, 0];
    let divInc = 360 / radialdivision;
    let triSplit = divInc / 2;
    let alpha = 0;

    let pt2Rad = [r*Math.cos(radians(alpha)), -0.5, r*Math.sin(radians(alpha))];
    while(alpha <= 360){
        //handle radial
        alpha += triSplit;
        let pt1Rad = [r*Math.cos(radians(alpha)), -0.5, r*Math.sin(radians(alpha))];
        addTriangle(ptCenterRad[0], ptCenterRad[1], ptCenterRad[2], pt2Rad[0], pt2Rad[1], pt2Rad[2], pt1Rad[0], pt1Rad[1], pt1Rad[2]);

        //handle sides
        let p1 = pt1Rad;
        let p2 = pt2Rad;
        //setup sides
        let distP1 = calcDist(p1, ptApex);
        let distP2 = calcDist(p2, ptApex);
        let xIncP1 = distP1[0] / heightdivision;
        let yIncP1 = distP1[1] / heightdivision;
        let zIncP1 = distP1[2] / heightdivision;
        let xIncP2 = distP2[0] / heightdivision;
        let yIncP2 = distP2[1] / heightdivision;
        let zIncP2 = distP2[2] / heightdivision;
        
        let runCnt = 0;
        while(runCnt < heightdivision){
            runCnt++;
            let p3 = [p1[0] + xIncP1, p1[1] + yIncP1, p1[2] + zIncP1];
            let p4 = [p2[0] + xIncP2, p2[1] + yIncP2, p2[2] + zIncP2];
            if(runCnt == heightdivision){
                addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], ptApex[0], ptApex[1], ptApex[2]);
            }
            else{
                addTriangle(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]);
                addTriangle(p2[0], p2[1], p2[2], p4[0], p4[1], p4[2], p3[0], p3[1], p3[2]);
            }
            p1 = p3;
            p2 = p4;
        }
        pt2Rad = pt1Rad;
    }
}

//calculates distance between two points for each coordinate individually
function calcDist(p1, p2){
    let xDist = Math.sqrt(Math.pow((p1[0] - p2[0]), 2));
    if(p2[0] < p1[0]){
        xDist = xDist - 2*xDist;
    }
    let yDist = Math.sqrt(Math.pow((p1[1] - p2[1]), 2));
    if(p2[1] < p1[1]){
        yDist = yDist - 2*yDist;
    }
    let zDist = Math.sqrt(Math.pow((p1[2] - p2[2]), 2));
    if(p2[2] < p1[2]){
        zDist = zDist - 2*zDist;
    }
    return [xDist, yDist, zDist];
}
    
//
// fill in code that creates the triangles for a sphere with diameter 1
// (centered at the origin) with number of slides (longitude) given by
// slices and the number of stacks (lattitude) given by stacks.
// For this function, you will implement the tessellation method based
// on spherical coordinates as described in the video (as opposed to the
//recursive subdivision method).
//
function makeSphere (subdivisions) {
    //make icosahedron
    let a = 2 / (1 + Math.sqrt(5));
    let v0 = normalizePoint([0, a, -1]);
    let v1 = normalizePoint([0-a, 1, 0]);
    let v2 = normalizePoint([a, 1, 0]);
    let v3 = normalizePoint([0, a, 1]);
    let v4 = normalizePoint([-1, 0, a]);
    let v5 = normalizePoint([0, 0-a, 1]);
    let v6 = normalizePoint([1, 0, a]);
    let v7 = normalizePoint([1, 0, 0-a]);
    let v8 = normalizePoint([0, 0-a, -1]);
    let v9 = normalizePoint([-1, 0, 0-a]);
    let v10 = normalizePoint([0-a, -1, 0]);
    let v11 = normalizePoint([a, -1, 0]);

    let tri0 = [v0, v1, v2];
    splitTriangle(tri0, subdivisions);
    let tri1 = [v3, v2, v1];
    splitTriangle(tri1, subdivisions);
    let tri2 = [v3, v4, v5];
    splitTriangle(tri2, subdivisions);
    let tri3 = [v3, v5, v6];
    splitTriangle(tri3, subdivisions);
    let tri4 = [v0, v7, v8];
    splitTriangle(tri4, subdivisions);
    let tri5 = [v0, v8, v9];
    splitTriangle(tri5, subdivisions);
    let tri6 = [v5, v10, v11];
    splitTriangle(tri6, subdivisions);
    let tri7 = [v8, v11, v10];
    splitTriangle(tri7, subdivisions);
    let tri8 = [v1, v9, v4];
    splitTriangle(tri8, subdivisions);
    let tri9 = [v10, v4, v9];
    splitTriangle(tri9, subdivisions);
    let tri10 = [v2, v6, v7];
    splitTriangle(tri10, subdivisions);
    let tri11 = [v11, v7, v6];
    splitTriangle(tri11, subdivisions);
    let tri12 = [v3, v1, v4];
    splitTriangle(tri12, subdivisions);
    let tri13 = [v3, v6, v2];
    splitTriangle(tri13, subdivisions);
    let tri14 = [v0, v9, v1];
    splitTriangle(tri14, subdivisions);
    let tri15 = [v0, v2, v7];
    splitTriangle(tri15, subdivisions);
    let tri16 = [v8, v10, v9];
    splitTriangle(tri16, subdivisions);
    let tri17 = [v8, v7, v11];
    splitTriangle(tri17, subdivisions);
    let tri18 = [v5, v4, v10];
    splitTriangle(tri18, subdivisions);
    let tri19 = [v5, v11, v6];
    splitTriangle(tri19, subdivisions);
}

function splitTriangle(tri, divs){
    if(divs < 2){
        addTriangle(tri[0][0], tri[0][1], tri[0][2], tri[1][0], tri[1][1], tri[1][2], tri[2][0], tri[2][1], tri[2][2]);
    }
    else{
        let mid1 = normalizePoint(calcMidpoint(tri[0], tri[1]));
        let mid2 = normalizePoint(calcMidpoint(tri[0], tri[2]));
        let mid3 = normalizePoint(calcMidpoint(tri[1], tri[2]));
        splitTriangle([mid1, mid3, mid2], divs-1);
        splitTriangle([tri[0], mid1, mid2], divs-1);
        splitTriangle([mid1, tri[1], mid3], divs-1);
        splitTriangle([mid2, mid3, tri[2]], divs-1);
    }
}

function normalizePoint(pt){
    let norm = Math.sqrt(Math.pow(pt[0], 2) + Math.pow(pt[1], 2) + Math.pow(pt[2], 2));
    return [pt[0]/norm, pt[1]/norm, pt[2]/norm];
}

////////////////////////////////////////////////////////////////////
//
//  Do not edit below this line
//
///////////////////////////////////////////////////////////////////

function radians(degrees)
{
  var pi = Math.PI;
  return degrees * (pi/180);
}

function addTriangle (x0,y0,z0,x1,y1,z1,x2,y2,z2) {

    
    var nverts = points.length / 4;
    
    // push first vertex
    points.push(x0);  bary.push (1.0);
    points.push(y0);  bary.push (0.0);
    points.push(z0);  bary.push (0.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++;
    
    // push second vertex
    points.push(x1); bary.push (0.0);
    points.push(y1); bary.push (1.0);
    points.push(z1); bary.push (0.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++
    
    // push third vertex
    points.push(x2); bary.push (0.0);
    points.push(y2); bary.push (0.0);
    points.push(z2); bary.push (1.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++;
}

