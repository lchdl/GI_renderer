#include "basedefs.h"

int getNumProcessors() {
    long nprocs = -1;
#ifdef _WIN32
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    nprocs = info.dwNumberOfProcessors;
#else /* LINUX */
    nprocs = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    return nprocs < 1 ? 1 : nprocs;
}

bool chInStr(const char ch, const char* s) {
    for (int i = 0; s[i] != '\0'; i++) {
        if (ch == s[i]) return true;
    }
    return false;
}
void strRemCh(const char * src, char * buf, int bufLen, char ch)
{
    memset(buf, 0, bufLen);
    int _s = 0, _t = 0;
    while (_s < strlen(src)) {
        if (src[_s] == ch) {
            _s++;
            continue;
        }
        if (_t >= bufLen-1)
            break;
        buf[_t++] = src[_s++];
    }
}
bool txtGetWord(FILE * fp, char * buf, int bufLen, const char* wordDelim, const char commentChar)
{
    if (feof(fp)) return false;
    int i = 0;
    int write = 0, end = 0;
    while (!feof(fp) && !end) {
        int c = fgetc(fp);

        if (chInStr(char(c), wordDelim) || c == -1) {
            if (c == commentChar) { /* skip comment */
                while ((c = fgetc(fp)) != EOF)
                    if (c == '\n')
                        break;
            }
            else if (write == 1) { /* a word is finished because it meets a delimiter */
                write = 0;
                end = 1;
            }
        }
        else { /* word still not complete, write to buffer */
            write = 1;
        }

        if (write && i < bufLen - 1) {
            buf[i] = char(c);
            i++;
        }
    }
    buf[i] = '\0'; /* reserve last character as '\0' */
    if (i == 0) return false; /* if parsed word is empty string return false */
    else return true;
}
bool strGetWord(char * src, char * buf, int bufLen, const char * wordDelim, const char commentChar)
{
    if (!src) return false;
    int i = 0;
    int s = 0, l = int(strlen(src));
    int write = 0, end = 0;
    while (s < l && !end) {
        int c = src[s++];

        if (chInStr(char(c), wordDelim)) {
            if (c == commentChar) { /* skip comment */
                while ((c = src[s++]) && (s < l))
                    if (c == '\n')
                        break;
            }
            else if (write == 1) { /* a word is finished because it meets a delimiter */
                write = 0;
                end = 1;
            }
        }
        else { /* word still not complete, write to buffer */
            write = 1;
        }

        if (write && i < bufLen - 1) {
            buf[i] = char(c);
            i++;
        }
    }
    buf[i] = '\0'; /* reserve last character as '\0' */
    /* offset source string by s characters */
    int j = 0;
    for (; s + j < l; j++) {
        src[j] = src[s + j];
    }
    src[j] = '\0';
    /* if parsed word is empty string return false */
    if (i == 0) return false;
    else return true;
}
bool txtGetLine(FILE * fp, char * buf, int bufLen)
{
    if (feof(fp)) return false;
    int i = 0;
    while (!feof(fp)) {
        int ch = fgetc(fp);
        if (ch == '\n') break;
        if (i >= bufLen - 1) break;
        buf[i] = char(ch);
        i++;
    }
    buf[i] = '\0';
    return true;
}
bool txtGetReal(FILE * fp, REAL * v)
{
    char* errptr = NULL, buf[32];
    if (!txtGetWord(fp, buf, 32, " #\n\r\t/", '#')) return false;
    *v = (REAL)(strtod(buf, &errptr));
    if (*errptr != '\0') return false;
    else return true;
}
bool txtGetInt(FILE * fp, int * v)
{
    char* errptr = NULL, buf[32];
    if (!txtGetWord(fp, buf, 32, " #\n\r\t/", '#')) return false;
    *v = (int)(strtol(buf, &errptr, 10));
    if (*errptr != '\0') return false;
    else return true;
}
bool txtGetVec2(FILE * fp, VEC2 * v)
{
    if (!txtGetReal(fp, &(v->x))) return false;
    if (!txtGetReal(fp, &(v->y))) return false;
    return true;
}
bool txtGetVec3(FILE * fp, VEC3 * v)
{
    if (!txtGetReal(fp, &(v->x))) return false;
    if (!txtGetReal(fp, &(v->y))) return false;
    if (!txtGetReal(fp, &(v->z))) return false;
    return true;
}
bool txtGetInt2(FILE * fp, INT2 * v)
{
    if (!txtGetInt(fp, &(v->x))) return false;
    if (!txtGetInt(fp, &(v->y))) return false;
    return true;
}
bool txtGetInt3(FILE * fp, INT3 * v)
{
    if (!txtGetInt(fp, &(v->x))) return false;
    if (!txtGetInt(fp, &(v->y))) return false;
    if (!txtGetInt(fp, &(v->z))) return false;
    return true;
}

bool loadTextFile(char * file, char * buf, int bufLen) {
    FILE* fp = NULL;
    fp = fopen(file, "rb");
    if (fp == NULL) {
        return false;
    }
    else {
        /* get file size */
        fseek(fp, 0, SEEK_END);
        int flen = int(ftell(fp));
        fseek(fp, 0, SEEK_SET);
        /* load file into mem */
        int readlen = flen > (bufLen - 1) ? (bufLen - 1) : flen;
        bool success = true;
        if (buf) {
            memset(buf, 0, bufLen);
            int len = int(fread(buf, 1, readlen, fp));
            if (len != readlen) { success = false; }
            else { success = true; }
            /* add '\0' to the end of the string */
            buf[readlen + 1] = '\0';
        }
        else { success = false; }
        fclose(fp);
        return success;
    }
}
Array<BYTE> loadAsByteArray(const char * file)
{
    Array<BYTE> byteArray;
    FILE* fp;
    if ((fp = fopen(file, "rb")) == NULL)
        return byteArray;
    BYTE byte;
    while (true) {
        if (fread(&byte, 1, 1, fp) != 1)
            break;
        else
            byteArray.append(byte);
    }
    fclose(fp);
    return byteArray;
}
bool saveByteArray(const Array<BYTE>& byteArray, const char * file)
{
    FILE* fp = NULL;
    if ((fp = fopen(file, "wb")) == NULL)
        return false;
    for (int i = 0; i < byteArray.size(); i++) {
        if (fwrite(&(byteArray[i]), 1, 1, fp) != 1){
            fclose(fp);
            return false;
        }
    }
    fclose(fp);
    return true;
}
void printVec3(const char* prefix, VEC3 v)
{
    if (prefix) {
        printf("%s(%.3lf, %.3lf, %.3lf)\n", prefix ,double(v.x), double(v.y), double(v.z));
    }
    else {
        printf("(%.3lf, %.3lf, %.3lf)\n", double(v.x), double(v.y), double(v.z));
    }
}
void printVec4(const char* prefix, VEC4 v)
{
    if (prefix) {
        printf("%s(%.3lf, %.3lf, %.3lf, %.3lf)\n", prefix, double(v.s), double(v.x), double(v.y), double(v.z));
    }
    else {
        printf("(%.3lf, %.3lf, %.3lf, %.3lf)\n", double(v.s), double(v.x), double(v.y), double(v.z));
    }
}
void printQuaternion(const char * prefix, QUAT q)
{
    if (prefix != NULL) {
        printf("%s", prefix);
    }
    printf("%.2f + %.2fi + %.2fj + %.2fk (%.2f, [%.2f, %.2f, %.2f])\n", 
        q.s, q.x, q.y, q.z, q.s, q.x, q.y, q.z);
}
void printByteArray(const Array<BYTE>& byteArray)
{
    for (int i = 0; i < byteArray.size(); i++)
        printf("%c", char(byteArray[i]));
}
void printMat3x3(const char * prefix, MAT3x3 m)
{
    if (prefix != NULL) {
        printf("%s\n", prefix);
    }
    printf("%.2f %.2f %.2f\n", m.xx, m.xy, m.xz);
    printf("%.2f %.2f %.2f\n", m.yx, m.yy, m.yz);
    printf("%.2f %.2f %.2f\n", m.zx, m.zy, m.zz);
}
void printMat4x4(const char * prefix, MAT4x4 m)
{
    if (prefix != NULL) {
        printf("%s\n", prefix);
    }
    printf("%.2f %.2f %.2f %.2f\n", m.xx, m.xy, m.xz, m.xs);
    printf("%.2f %.2f %.2f %.2f\n", m.yx, m.yy, m.yz, m.ys);
    printf("%.2f %.2f %.2f %.2f\n", m.zx, m.zy, m.zz, m.zs);
    printf("%.2f %.2f %.2f %.2f\n", m.sx, m.sy, m.sz, m.ss);
}

String::String()
{
    append('\0');
}

String::String(const String & that)
{
    append(that);
}

String::String(const char * s)
{
    for (int i = 0; i < int(strlen(s)); i++)
        append(s[i]);
    append('\0');
}

String String::operator=(const String & that)
{
    if (this == &that)
        return (*this);
    clear();
    append(that);

    return (*this);
}

String String::operator=(const char * s)
{
    clear();

    if (s != NULL)
        for (int i = 0; i<int(strlen(s)); i++)
            append(s[i]);
    append('\0');

    return (*this);
}

String::~String()
{
}

bool String::operator==(const char * s)
{
    int l = int(strlen(s));
    if (length() != l) return false;
    for (int i = 0; i < l; i++) {
        if (baseptr[i] != s[i])
            return false;
    }
    return true;
}

bool String::operator==(const String & that)
{
    if (this->length() != that.length()) return false;
    for (int i = 0; i < that.length(); i++) {
        if (baseptr[i] != that.baseptr[i])
            return false;
    }
    return true;
}

int String::length() const
{
    return this->Ne - 1;
}

const char * String::c_str() const
{
    return baseptr;
}

String::operator const char* () const
{
    return c_str();
}

String String::operator+(const String & s) const
{
    if (this->length() == 0)
        return s;
    if (s.length() == 0)
        return (*this);
    
    String newString;
    newString[0] = this->at(0); /* remove '\0' */
    for (int i = 1; i < length(); i++)
        newString.append(this->at(i));
    for (int i = 0; i < s.length(); i++)
        newString.append(s[i]);
    newString.append('\0');

    return newString;
}

String String::operator+=(const String & s)
{
    if (s.length() == 0)
        return (*this);

    this->at(length()) = s[0];
    for (int i = 1; i < s.length(); i++)
        append(s[i]);
    this->append('\0');

    return (*this);
}
