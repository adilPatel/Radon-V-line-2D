/* SLLineIntegralApproximator
 * The line integral approximator previously coded in approximate_radon,
 * but a super legera version built in C. This function will approximate a
 * line integral in a square voxel domain.
 */

#include "mex.h"
#include "matrix.h"
#include <math.h>

#define NO_INTERSECT -2.0

/*
 * Structure Vector3: An ordinary 3-component vector.
 */
typedef struct Vector3 {
    double x;
    double y;
    double z;
} Vector3;

/*
 * Enumeration AxisIndex: Assigns each axis to a number.
 */
typedef enum AxisIndex {
    x_axis = 0,
    y_axis = 1,
    z_axis = 2
} AxisIndex;

/*
 * Structure AxisTime: Tags a time value with its associated axis.
 */
typedef struct AxisTime {
    AxisIndex axis;
    double t;
} AxisTime;

/*
 * Structure TimeValues: Contains the sorted time intersections, tagged
 *                       their corresponding axes.
 */
typedef struct TimeValues {
    AxisTime first;
    AxisTime second;
    AxisTime third;
} TimeValues;

TimeValues closestIntersection(Vector3 r_0, Vector3 v, Vector3 planes);
Vector3 calculateExcess(Vector3 r_0, Vector3 v, const mxArray *grid_values, int grid_width, double &output, TimeValues &intersections);
// The necessary vector function declarations.
double vectorNorm(Vector3 vec);
Vector3 makeVector3(double x, double y, double z);
Vector3 normaliseVector(Vector3 vec);
Vector3 scalarMultiplyVector(Vector3 vec, double scalar);
Vector3 addVectors(Vector3 vec1, Vector3 vec2);
Vector3 subtractVectors(Vector3 vec1, Vector3 vec2);
Vector3 initVector(double *vector);

int contiguousArrayIndex(const mxArray *array, Vector3 vec, int grid_width);
void validateArguments(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]);

/*
 * Function mexFunction: Acts as the main entry point of the function.
 * Parameter nlhs: The number of output arguments.
 * Parameter plhs: The output data of the function.
 * Parameter nrhs: The number of input arguments.
 * Parameter prhs: The input data of the function.
 */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    // Validate everything before anything.
    validateArguments(nlhs, plhs, nrhs, prhs);
    
    // The values of the discretised 4D function.
    const mxArray *grid_values = prhs[0];
    
    const double *grid_doubles = mxGetPr(grid_values);
    
    // The distance from the centre of the grid to the edge of the domain.
    int grid_width  = (int)mxGetScalar(prhs[1]);
    
    // The initial point for the ray vector.
    Vector3 r_0 = initVector(mxGetPr(prhs[2]));
    
    // The final point for the ray vector.
    Vector3 r_1 = initVector(mxGetPr(prhs[3]));
    
    // Self-explanatory. I'd be very worried if you didn't understand the
    // name of this variable.
    double output = 0.0;
    
    // First, we retrieve the velocity vector and normalise it...
    const Vector3 v = normaliseVector(subtractVectors(r_1, r_0));
    
    // Find the excess in the beginning, plus calculate the starting point.
    
    TimeValues start_int;
    Vector3 start = calculateExcess(r_0, v, grid_values, grid_width, output, start_int);
    
    // We also need to do something similar with r_1. An easy way to tackle
    // this is to simply execute the above line with the v direction
    // flipped, and with r_1 behaving like r_0...
    TimeValues end_int;
    const Vector3 end = calculateExcess(r_1, v, grid_values, grid_width, output, end_int);
    
    // Nice and clean when putting it in a function. I just wish there was
    // code folding in MATLAB. Maybe I'll link them to separate files one
    // day.
    
    // Now for the main show. Following this is a loop where we go through
    // all the voxels along the ray and add them up. Following the previous
    // section, we want to find the planes that are in the direction of the
    // v vector. We will then compute the distance that the ray has in the
    // voxel and multiply it with the voxel value. Meanwhile, we accumulate
    // everything. This is the climax of the function... we'll just be
    // cleaning up after this...
    
    
    double reference_length = vectorNorm(subtractVectors(r_1, r_0));
    
    Vector3 current_pos = start;
    
    int x_add_factor = (v.x > 0) ? 1 : 0;
    int y_add_factor = (v.y > 0) ? 1 : 0;
    int z_add_factor = (v.z > 0) ? 1 : 0;
    
    int x_test = (int)floor(start.x) + x_add_factor;
    int y_test = (int)floor(start.y) + y_add_factor;
    int z_test = (int)floor(start.z) + z_add_factor;
    
    x_test += (fmod(r_0.x, 1.0) == 0.0) ? x_add_factor - 1 : 0;
    y_test += (fmod(r_0.y, 1.0) == 0.0) ? y_add_factor - 1 : 0;
    z_test += (fmod(r_0.z, 1.0) == 0.0) ? z_add_factor - 1 : 0;
    
    while (true) {
             
      Vector3 test_planes = makeVector3((double)x_test, (double)y_test, (double)z_test);
        
        // Test the intersections...
        TimeValues intersections = closestIntersection(r_0, v, test_planes);
        
        // Calculate the next point...
        double t = intersections.first.t;
        Vector3 next_pos = addVectors(r_0, scalarMultiplyVector(v, t));
        
        // These lines will account for the moments when the ray hits a 
        // corner. It might be overengineered, but it's robust (hopefully).
        bool overlap12  = intersections.first.t == intersections.second.t;
        bool overlap23  = intersections.second.t == intersections.third.t;
        bool overlap123 = overlap12 && overlap23;
        
        AxisIndex first_ax  = intersections.first.axis;
        AxisIndex second_ax = intersections.second.axis;
        AxisIndex third_ax  = intersections.third.axis;
        
        bool x_cond = (first_ax == x_axis) || ((second_ax == x_axis) && overlap12) || overlap123;
        bool y_cond = (first_ax == y_axis) || ((second_ax == y_axis) && overlap12) || overlap123;
        bool z_cond = (first_ax == z_axis) || ((second_ax == z_axis) && overlap12) || overlap123;
        
        x_test += x_cond ? x_add_factor - 1 : 0;
        y_test += y_cond ? y_add_factor - 1 : 0;
        z_test += z_cond ? z_add_factor - 1 : 0;

        
        // First, we need to check if we have to terminate the loop. We do this
        // by checking if the total length (including the next intersection)
        // exceeds the original length...
        double total_length = vectorNorm(subtractVectors(next_pos, r_0));
        
        if (total_length > reference_length) {
            break;
        }
        
        // Perform the approximation...
        double norm = vectorNorm(subtractVectors(next_pos, current_pos));
        int i = contiguousArrayIndex(grid_values, current_pos, grid_width);
        double voxel_value = grid_doubles[i];
                
        current_pos = next_pos;
        
        output += voxel_value * norm;
        
        x_test += x_add_factor;
        y_test += y_add_factor;
        z_test += z_add_factor;
        
    }
    
    plhs[0] = mxCreateDoubleScalar(output);
    
}

/*
 * Function calculateExcess: If the ray starts or ends in the middle of a
 *                           voxel, this moves the origin of the ray to a 
 *                           plane and adds the contribution from this 
 *                           shift.
 * Parameter r_0: The start of the ray.
 * Parameter v: The velocity vector of the ray.
 * Parameter grid_values: The discrete function grid.
 * Parameter grid_width: The width of the domain from the axis to its 
 *                       origin.
 * Parameter output: The output of the  mexFunction; the result of the
 *                   integral 
 * Returns: A structure containing sorted intersection times and their 
 *          corresponding axes.
 */
Vector3 calculateExcess(Vector3 r_0, Vector3 v, const mxArray *grid_values, int grid_width, double &output, TimeValues &intersections) {
    
    // Calculate the nearest planes to test for intersection...
    double x_test, y_test, z_test;
    
    x_test = (v.x > 0) ? floor(r_0.x) + 1.0 : floor(r_0.x);
    y_test = (v.y > 0) ? floor(r_0.y) + 1.0 : floor(r_0.y);
    z_test = (v.z > 0) ? floor(r_0.z) + 1.0 : floor(r_0.z);
    // Can't do these shortcuts with MATLAB, can you?
    
    Vector3 test_planes = makeVector3(x_test, y_test, z_test);
    
    // Below are the the sorted time values with their associated axis
    // indices for easy reference.
    intersections = closestIntersection(r_0, v, test_planes);
    
    // Btw all this crap with the sorted and tagged structures are used to
    // sort out all the pathological cases like r_0 starting in a junction
    // between 8 voxels. Bloody pathological cases always spoil the fun >:(
    
    // We now add in the partial part (if there is any...)
    
    const double *grid_doubles = mxGetPr(grid_values);
    
    if (intersections.first.t > 0.0) {
        
        double t_0 = intersections.first.t;
        Vector3 nearest_r0 = addVectors(r_0, scalarMultiplyVector(v, t_0));
        double ds = vectorNorm(subtractVectors(nearest_r0, r_0));
        int i = contiguousArrayIndex(grid_values, r_0, grid_width);
        
        output += grid_doubles[i] * ds;
        return nearest_r0;
        
    }
    
   return r_0;
    
}

/*
 * Function closestIntersections: Given a plane, this will find the nearest
 *                                intersection time and axis.
 * Parameter r_0: The start of the ray.
 * Parameter v: The velocity vector of the ray.
 * Parameter planes: A vector containing the axes you want to test the 
 *                   intersections of.
 * Returns: A structure containing sorted intersection times and their 
 *          corresponding axes.
 */
TimeValues closestIntersection(Vector3 r_0, Vector3 v, Vector3 planes) {
    
    // We mark the time as NO_INTERSECT if the v component never intersects
    // the plane...
    double t_x = (v.x != 0.0) ? (planes.x - r_0.x) / v.x : NO_INTERSECT;
    double t_y = (v.y != 0.0) ? (planes.y - r_0.y) / v.y : NO_INTERSECT;
    double t_z = (v.z != 0.0) ? (planes.z - r_0.z) / v.z : NO_INTERSECT;
        
    // These are just compact ways of marking intersections. They'll be 
    // useful later.
    bool x_int = t_x != NO_INTERSECT;
    bool y_int = t_y != NO_INTERSECT;
    bool z_int = t_z != NO_INTERSECT;
    
    // Tagging the axes with their times here...
    AxisTime a_x, a_y, a_z;
    a_x.axis = x_axis;
    a_y.axis = y_axis;
    a_z.axis = z_axis;
    a_x.t = t_x;
    a_y.t = t_y;
    a_z.t = t_z;
    
    // Calculate and sort the minima. At this point, I hope you're ready
    // for one of the most confusing lines of code you've ever seen.
    // This part calculates and sorts the time value from minimum to
    // maximum. The dirty uses of the ternary operators work as follows:
    // For the first line, calculate the minimum between X and Y. It's
    // straightforward if both the X and Y components have finite
    // plane intersections (i.e. non-zero v components). We just return
    // the X if the X time <= Y time. You can see it here. But what if X
    // has no intersections? The flag NO_INTERSECTIONS is negative so we
    // have to account for that. Thus, we throw it into another ternary
    // operator that checks first if it has finite intersections. If it
    // doesn't, we return the Y axis. What if Y has no intersections, but
    // X does? Again, in the nested ternary operator, we check if X has
    // finite intersections. But what if none have finite intersections?
    // In this case, it doesn't matter at all then; we discard both. But
    // what if Y <= X while both have finite intersections? Well, the next
    // line should answer your question.
    AxisTime min_xy = ((a_x.t <= a_y.t) && x_int && y_int) ? a_x : (x_int ? a_x : a_y);
    min_xy = ((a_y.t <= a_x.t) && x_int && y_int) ? a_y : min_xy;
    bool xy_int = min_xy.t != NO_INTERSECT;
    AxisTime first  = ((a_z.t <= min_xy.t) && xy_int && z_int) ? a_z : (z_int ? a_z : min_xy);
    first = ((min_xy.t <= a_z.t) && xy_int && z_int) ? min_xy : first;
    AxisTime second = (a_z.t == first.t) ? min_xy : a_z;
    AxisTime third  = (a_x.t == min_xy.t) ? a_y : a_x;
        
    // Store it all in a structure that we pack up and send away!
    TimeValues sorted_values;
    sorted_values.first  = first;
    sorted_values.second = second;
    sorted_values.third  = third;
    
    return sorted_values;
    
}

Vector3 makeVector3(double x, double y, double z) {
    Vector3 output;
    output.x = x;
    output.y = y;
    output.z = z;
    return output;
}

double vectorNorm(Vector3 vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

Vector3 normaliseVector(Vector3 vec) {
    double normFactor = 1.0 / vectorNorm(vec);
    return scalarMultiplyVector(vec, normFactor);
}

Vector3 scalarMultiplyVector(Vector3 vec, double scalar) {
    Vector3 output;
    output.x = vec.x * scalar;
    output.y = vec.y * scalar;
    output.z = vec.z * scalar;
    return output;
}

Vector3 subtractVectors(Vector3 vec1, Vector3 vec2) {
    Vector3 output;
    output.x = vec1.x - vec2.x;
    output.y = vec1.y - vec2.y;
    output.z = vec1.z - vec2.z;
    return output;
}

Vector3 addVectors(Vector3 vec1, Vector3 vec2) {
    Vector3 output;
    output.x = vec1.x + vec2.x;
    output.y = vec1.y + vec2.y;
    output.z = vec1.z + vec2.z;
    return output;
}

Vector3 initVector(double *vector) {
    Vector3 output;
    output.x = vector[0];
    output.y = vector[1];
    output.z = vector[2];
    return output;
}

/*
 * Function contiguousArrayIndex: Converts a 3D array index into a 1D
 *                                contiguous memory index.
 * Parameter array: The 3D array.
 * Parameter vec: The position vector you want to find the indices of.
 * Parameter grid_width: The distance from the origin to the edge of the
 *                       domain.
 * Returns: The 1D contiguous memory index.
 */
int contiguousArrayIndex(const mxArray *array, Vector3 vec, int grid_width) {
    
    int x_index = (int)vec.x + grid_width;
    int y_index = (int)vec.y + grid_width;
    int z_index = (int)vec.z + grid_width;
        
    mwIndex *subs = (mwIndex *)mxCalloc(3, sizeof(mwIndex));
    subs[0] = (mwIndex)x_index;
    subs[1] = (mwIndex)y_index;
    subs[2] = (mwIndex)z_index;
    
    return (int)mxCalcSingleSubscript(array, 3, subs);
    
    
}

void validateArguments(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    
    // First, validate the number of input arguments, as well as their
    // dimensions...
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:SLLineIntegralApproximator:nrhs",
                "Four inputs required.");
    }
    
    // Now the number of output arguments...
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:SLLineIntegralApproximator:nlhs",
                "One output required.");
    }
    
    // After checking the number of arguments, we will check their
    // respective dimensions. I will show the actual arguments later on...
    if (mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:SLLineIntegralApproximator:notScalar",
                "Input 2 must be a scalar.");
    }
    
    
    if (mxGetN(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("MyToolbox:SLLineIntegralApproximator:notScalar",
                "Input 3 must be a vector.");
    }
    
    if (mxGetN(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("MyToolbox:SLLineIntegralApproximator:notScalar",
                "Input 4 must be a vector.");
    }
    
}