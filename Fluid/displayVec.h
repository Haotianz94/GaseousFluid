#ifndef _DISPLAYVEC_H_
#define _DISPLAYVEC_H_

#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>

class DisplayVec
{
private:
	int L;
	cv::Mat input;
	const char *inputName;
	int wide;
	int height;
	int size;
	int *pixel;

	float LIC_integral(float a, float b);

public:
	DisplayVec(int l) :L(l) {}
	DisplayVec(int l, int w, int h);
	~DisplayVec();
	void createInputTexture();
	void loadInputTexture();
	float getOutputTextureDDA(int x, int y, float *Vx, float *Vy);
	float getOutputTextureLIC(int x, int y, float *Vx, float *Vy);
	
	void test(const char *, bool DDA);
	void testDDA(int x, int y, float k);
};

#endif