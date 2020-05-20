#if !defined(VX3_MATHTREE_H)
#define VX3_MATHTREE_H

enum VX3_MathTreeOperator {
    mtEND,
    mtCONST,
    mtE,  // number e
    mtPI, // number pi
    mtVAR,
    mtADD,
    mtSUB,
    mtMUL,
    mtDIV,
    mtPOW,  // power
    mtSQRT, // sqrt
    mtSIN,
    mtCOS,
    mtTAN,
    mtATAN,
    mtLOG, // log_e
    mtINT, // round to nearest integer. e.g. 0.9 --> 1.0
    mtABS,
    mtNOT,
    mtGREATERTHAN,
    mtLESSTHAN,
    mtAND,
    mtOR,
    mtNORMALCDF, // normal CDF function
};
struct VX3_MathTreeToken {
    VX3_MathTreeOperator op = mtEND;
    double value = 0.0;
    void set(VX3_MathTreeOperator inOp, double inValue = 0.0) {
        op = inOp;
        value = inValue;
    }
};
struct VX3_MathTree {
    static bool validate(VX3_MathTreeToken *buff) {
        try {
            eval(1, 1, 1, 1, 1, 1, 1, 1, 1, buff);
        } catch (...) {
            return false;
        }
        return true;
    }
    /* Standard implementation is CUDA Math API
    https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html
    */
    __device__ __host__ static double eval(double x, double y, double z, double hit, double t, double angle, double closeness,
                                           int numClosePairs, int num_voxel, VX3_MathTreeToken *buff) {
        double values[1024];
        int values_cursor = 0;
        int process_cursor = 0;
        for (int i = 0; i < 1024; i++) {
            switch (buff[i].op) {
            case mtEND:
                return values[process_cursor];
            case mtCONST:
                values[values_cursor] = buff[i].value;
                break;
            case mtE:
                values[values_cursor] = 2.71828182845904523536;
                break;
            case mtPI:
                values[values_cursor] = 3.14159265358979323846;
                break;
            case mtVAR:
                if (buff[i].value < 0.5) {
                    values[values_cursor] = x;
                } else if (buff[i].value < 1.5) {
                    values[values_cursor] = y;
                } else if (buff[i].value < 2.5) {
                    values[values_cursor] = z;
                } else if (buff[i].value < 3.5) {
                    values[values_cursor] = hit;
                } else if (buff[i].value < 4.5) {
                    values[values_cursor] = t;
                } else if (buff[i].value < 5.5) {
                    values[values_cursor] = angle;
                } else if (buff[i].value < 6.5) {
                    values[values_cursor] = closeness;
                } else if (buff[i].value < 7.5) {
                    values[values_cursor] = numClosePairs;
                } else if (buff[i].value < 8.5) {
                    values[values_cursor] = num_voxel;
                } else {
                    // ERROR
                }
                break;
            case mtSIN:
                values[values_cursor] = sin(values[process_cursor]);
                process_cursor++;
                break;
            case mtCOS:
                values[values_cursor] = cos(values[process_cursor]);
                process_cursor++;
                break;
            case mtTAN:
                values[values_cursor] = tan(values[process_cursor]);
                process_cursor++;
                break;
            case mtATAN:
                values[values_cursor] = atan(values[process_cursor]);
                process_cursor++;
                break;
            case mtLOG:
                values[values_cursor] = log(values[process_cursor]);
                process_cursor++;
                break;
            case mtINT:
                values[values_cursor] = rint(values[process_cursor]);
                process_cursor++;
                break;
            case mtNORMALCDF:
                values[values_cursor] = normcdf(values[process_cursor]);
                process_cursor++;
                break;
            case mtADD:
                values[values_cursor] = values[process_cursor + 1] + values[process_cursor];
                process_cursor += 2;
                break;
            case mtSUB:
                values[values_cursor] = values[process_cursor + 1] - values[process_cursor];
                process_cursor += 2;
                break;
            case mtMUL:
                values[values_cursor] = values[process_cursor + 1] * values[process_cursor];
                process_cursor += 2;
                break;
            case mtDIV:
                values[values_cursor] = values[process_cursor + 1] / values[process_cursor];
                process_cursor += 2;
                break;
            case mtPOW:
                values[values_cursor] = pow(values[process_cursor + 1], values[process_cursor]);
                process_cursor += 2;
                break;
            case mtSQRT:
                values[values_cursor] = sqrt(values[process_cursor]);
                process_cursor++;
                break;
            case mtABS:
                values[values_cursor] = abs(values[process_cursor]);
                process_cursor++;
                break;
            case mtNOT:
                values[values_cursor] = !values[process_cursor];
                process_cursor++;
                break;
            case mtGREATERTHAN:
                values[values_cursor] = values[process_cursor + 1] > values[process_cursor];
                process_cursor += 2;
                break;
            case mtLESSTHAN:
                values[values_cursor] = values[process_cursor + 1] < values[process_cursor];
                process_cursor += 2;
                break;
            case mtAND:
                values[values_cursor] = values[process_cursor + 1] && values[process_cursor];
                process_cursor += 2;
                break;
            case mtOR:
                values[values_cursor] = values[process_cursor + 1] || values[process_cursor];
                process_cursor += 2;
                break;

            default:
#ifdef __CUDA_ARCH__
                printf("ERROR: not implemented.\n");
                return -1;
#else
                throw 0;
#endif
            }
            if (process_cursor > values_cursor) {
#ifdef __CUDA_ARCH__
                printf("ERROR: token sequence error: process_cursor > values_cursor.\n");
                return -1;
#else
                throw 1;
#endif
            }
            values_cursor++;
        }
#ifdef __CUDA_ARCH__
        printf("ERROR: VX3_MathTree Overflow.\n");
        return -1;
#else
        throw 2;
#endif
    }
};

#endif // VX3_MATHTREE_H
