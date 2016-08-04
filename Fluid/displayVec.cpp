#include "displayVec.h"
#include <cstdio>
#include <ctime>
#include <Eigen\Eigen>
#include <math.h>

#define IX(x, y) ( (x) + (y) * (wide + 2) )
#define BOUNDED(x, y) ( ((x) < 1 || (x) > wide+1 || (y) < 1 || (y) > height+1)? false : true)
#define PI 3.14159265

DisplayVec::DisplayVec(int l, int w, int h)
{
	L = l;
	wide = w;
	height = h;
	pixel = new int [(wide+2)*(height+2)]; 
	//for(int y = 0; y < height; y++)
	//	pixel[y] = new int [wide];
}

DisplayVec::~DisplayVec()
{
	//for(int y = 0; y < height; y++)
	//	delete [] pixel[y];
	delete[] pixel;
}

void DisplayVec::createInputTexture()
{
	srand(time(0));
	input = cv::Mat(height, wide, CV_8UC1);

	for(int y = 1; y <= height; y++)
		for(int x = 1; x <= wide; x++)
		{
			pixel[IX(x,y)] = rand()%256;
			input.at<uchar>(height-y, x-1) = pixel[IX(x,y)];
		}

	//box filter
	/*int step = 1;
	for(int y = 1; y <= height; y++)
		for(int x = 1; x <= wide; x++)
		{
			int sum = 0;
			for(int i = -step; i <= step; i++)
				for(int j = -step; j <= step; j++)
				{
					if(x+i >= 1 && x+i <= wide && y+j >= 1 && y+j <= height)
						sum += pixel[IX(x+i, y+j)];
					else
						sum += pixel[IX(x, y)];
				}
				pixel[IX(x, y)] = sum / ((2*step+1)*(2*step+1));
				input.at<uchar>(height-y, x-1) = pixel[IX(x,y)];
		}
	*/
	//cv::imshow("texture", input);
	//cv::waitKey(0);
}

void DisplayVec::loadInputTexture()
{
	input = cv::imread(inputName, 0);

	for(int y = 1; y <= height; y++)
		for(int x = 1; x <= wide; x++)
			pixel[IX(x, y)] = input.at<uchar>(height-y, x-1);
}

float DisplayVec::getOutputTextureDDA(int x, int y, float *Vx, float *Vy)
{
	float k = Vy[IX(x,y)] / Vx[IX(x,y)];
	float fk = fabsf(k);
	int dy = (k >= 0)? 1 : -1;
	int dx = dy;
	
	int num = 1;
	int integral = pixel[IX(x, y)];
	float e = -0.5;
	int xp = x, xn = x;
	int yp = y, yn = y;
	if(fk <= 1)
	{
		for(int i = 0; i < L; i++)
		{
			xp += 1;
			xn -= 1;
			e += fk;
			if(e >= 0)
			{
				yp += dy;
				yn -= dy;
				e -= 1;
			}
			if(xp <= wide && yp >= 1 && yp <= height)
			{
				integral += pixel[IX(xp, yp)];
				num ++;
			}
			if(xn >= 1 && yn >= 1 && yn <= height)
			{
				integral += pixel[IX(xn, yn)];
				num ++;
			}
		}
	}
	else
	{
		fk = 1.0 / fk;
		for(int i = 0; i < L; i++)
		{
			yp += 1;
			yn -= 1;
			e += fk;
			if(e >= 0)
			{
				xp += dx;
				xn -= dx;
				e -= 1;
			}
			if(yp <= height && xp >= 1 && xp <= wide)
			{
				integral += pixel[IX(xp, yp)];
				num ++;
			}
			if(yn >= 1 && xn >= 1 && xn <= wide)
			{
				integral += pixel[IX(xn, yn)];
				num ++;
			}
		}
	}

	return 1.0f * integral / num / 255;
}

void DisplayVec::testDDA(int x, int y, float k)
{
	float fk = fabsf(k);
	int dy = (k >= 0)? 1 : -1;
	int dx = dy;
	
	int num = 1;
	int integral = pixel[IX(x, y)];
	float e = -0.5;
	int xp = x, xn = x;
	int yp = y, yn = y;
	if(fk <= 1)
	{
		for(int i = 0; i < L; i++)
		{
			xp += 1;
			xn -= 1;
			e += fk;
			if(e >= 0)
			{
				yp += dy;
				yn -= dy;
				e -= 1;
			}
			if(xp <= wide && yp >= 1 && yp <= height)
			{
				integral += pixel[IX(xp, yp)];
				input.at<uchar>(height - yp, xp - 1) = 255;
				num ++;
			}
			if(xn >= 1 && yn >= 1 && yn <= height)
			{
				integral += pixel[IX(xn, yn)];
				input.at<uchar>(height - yn, xn - 1) = 255;
				num ++;
			}
		}
	}
	else
	{
		fk = 1.0 / fk;
		for(int i = 0; i < L; i++)
		{
			yp += 1;
			yn -= 1;
			e += fk;
			if(e >= 0)
			{
				xp += dx;
				xn -= dx;
				e -= 1;
			}
			if(yp <= height && xp >= 1 && xp <= wide)
			{
				integral += pixel[IX(xp, yp)];
				input.at<uchar>(height - yp, xp - 1) = 255;
				num ++;
			}
			if(yn >= 1 && xn >= 1 && xn <= wide)
			{
				integral += pixel[IX(xn, yn)];
				input.at<uchar>(height - yn, xn - 1) = 255;
				num ++;
			}
		}
	}
	
	cv::imshow("DDA", input);
	cv::waitKey(0);
}

float DisplayVec::getOutputTextureLIC(int x, int y, float *Vx, float *Vy)
{
	//int size = 10;
	//cv::Mat stream = cv::Mat(height*(size+1), wide*(size+1), CV_8UC3, cvScalar(0,0,0));

	cv::Point curGrid(x, y), lastGrid(-1, -1), nextGrid;
	cv::Point2f P(x+0.5, y+0.5), lastP(x+0.5, y+0.5);
	
	float ds, S = 0, s[4];
	float vx, vy;
	float MAX = 1e10;
	float eps = 0.001;
	float numerator = 0, denominator = 0;

	//forward 
	while(S <= L)
	{
  		vx = Vx[IX(curGrid.x, curGrid.y)];
		vy = Vy[IX(curGrid.x, curGrid.y)];
		
		if(vx == 0 && vy == 0)
			break;
		if(!BOUNDED(P.x, P.y))
			break;

		s[0] = (floorf(P.x) - P.x) / vx;
		if(s[0] < 0)
			s[0] = MAX;
		s[1] = (ceilf(P.x) - P.x) / vx;
		if(s[1] < 0)
			s[1] = MAX;
		s[2] = (ceilf(P.y) - P.y) / vy;
		if(s[2] < 0)
			s[2] = MAX;
		s[3] = (floorf(P.y) - P.y) / vy;
		if(s[3] < 0)
			s[3] = MAX;

		//find min 4
		ds = MAX;
		for(int i = 0; i < 4; i++)
			if(s[i] < ds)
				ds = s[i];
		ds += eps;
		P += cv::Point2f(vx, vy) * ds;
		//draw
		//cv::Point2f a = cv::Point2f((lastP.x-1)*size, (height+1-lastP.y)*size);
		//cv::Point2f b = cv::Point2f((P.x-1)*size, (height+1-P.y)*size);
		//cv::line(stream, a, b, cvScalar(255,255,255));
		//lastP = P;

		nextGrid = cv::Point(int(P.x), int(P.y));
		if(nextGrid == lastGrid)
			break;

		//calculate integral
		ds *= sqrtf(vx*vx + vy*vy);
		//Hanning window and ripple filter function
		//float h = LIC_integral(S, S+ds);
		//constant integral
		float h = ds;
		numerator += pixel[IX(curGrid.x, curGrid.y)] * h;
		denominator += h;
		
		S += ds;
		lastGrid = curGrid;
		curGrid = nextGrid;
	}

	curGrid = cv::Point(x, y);
	lastGrid = cv::Point(-1, -1);
	P = cv::Point2f(x+0.5, y+0.5);
	lastP = P;
	S = 0;
	while(S <= L)
	{
		vx = Vx[IX(curGrid.x, curGrid.y)];
		vy = Vy[IX(curGrid.x, curGrid.y)];
		
		if(vx == 0 && vy == 0)
			break;
		if(!BOUNDED(P.x, P.y))
			break;

		s[0] = (floorf(P.x) - P.x) / vx;
		if(s[0] < 0)
			s[0] = MAX;
		s[1] = (ceilf(P.x) - P.x) / vx;
		if(s[1] < 0)
			s[1] = MAX;
		s[2] = (ceilf(P.y) - P.y) / vy;
		if(s[2] < 0)
			s[2] = MAX;
		s[3] = (floorf(P.y) - P.y) / vy;
		if(s[3] < 0)
			s[3] = MAX;

		//find min 4
		ds = MAX;
		for(int i = 0; i < 4; i++)
			if(s[i] < ds)
				ds = s[i];
		ds += eps;
		P += cv::Point2f(vx, vy) * ds;
		//draw
		//cv::Point2f a = cv::Point2f((lastP.x-1)*size, (height+1-lastP.y)*size);
		//cv::Point2f b = cv::Point2f((P.x-1)*size, (height+1-P.y)*size);
		//cv::line(stream, a, b, cvScalar(255,255,255));
		//lastP = P;

		nextGrid = cv::Point(int(P.x), int(P.y));
		if(nextGrid == lastGrid)
			break;

		//calculate integral
		ds *= sqrtf(vx*vx + vy*vy); 
		//Hanning window and ripple filter function
		//float h = LIC_integral(S, S+ds);
		//constant integral
		float h = ds;
		numerator += pixel[IX(curGrid.x, curGrid.y)] * h;
		denominator += h;
		
		S += ds;
		lastGrid = curGrid;
		curGrid = nextGrid;
	}
	//stream.at<cv::Vec3b>((height-y+0.5)*size, (x-0.5)*size) = cv::Vec3b(0, 0, 255);
	//cv::imshow("stream line", stream);
	//cv::waitKey(0);

	return numerator / denominator / 255;
}

float DisplayVec::LIC_integral(float a, float b)
{
	float c = 2.9, d = 3, beita = PI/3;

	float res = 0;
	res += b - a + (sinf(b*c)-sinf(a*c)) / c;
	res += (sinf(b*d+beita) - sinf(a*d+beita)) / d;
	res += (sinf(b*(c-d)-beita) - sinf(a*(c-d)-beita)) / (2*(c-d));
	res += (sinf(b*(c+d)+beita) - sinf(a*(c+d)+beita)) / (2*(c+d));
	return res / 4;
}

/*void DisplayVec::drawArrow(int i, int j, float vx, float vy)
{
	float cx = (i+0.5) * GRIDSIZE;
	float cy = (j+0.5) * GRIDSIZE;
	float vl = sqrtf(vx*vx + vy*vy);
	//float seita = atanf(fabs(vy / vx));

	float max_v = max_vx > max_vy? max_vx : max_vy;
	if(max_v == 0)
		return;
	float dx = 0.5 * GRIDSIZE * vx / vl;
	float dy = 0.5 * GRIDSIZE * vy / vl;

	//float dx = 20 * GRIDSIZE * vx ;
	//float dy = 20 * GRIDSIZE * vy ;

	glBegin(GL_LINES);
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex2f(cx, cy);
	//glColor3f(0.0f, 0.0f, 1.0f);
	glVertex2f(cx + dx, cy + dy);
	
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex2f(cx, cy);
	glVertex2f(cx - vx/vl*GRIDSIZE*0.2, cy - vy/vl*GRIDSIZE*0.2);
	glEnd();
}*/

float DisplayVec::interpolate(int X, int Y, float *u)
{
	float x0 = 1.0 * X / gridSize - 0.5;
	float y0 = 1.0 * Y / gridSize - 0.5;

	int i0 = int(x0), i1 = i0 + 1;
	int j0 = int(y0), j1 = j0 + 1;
	float s1 = x0 - i0, s0 = 1 - s1;
	float t1 = y0 - j0, t0 = 1 - t1;

	return s0 * (t0*u[IX(i0,j0)] + t1*u[IX(i0,j1)]) +
		   s1 * (t0*u[IX(i1,j0)] + t1*u[IX(i1,j1)]);
}

void DisplayVec::test(const char *file, bool DDA)
{
	cv::Mat image = cv::imread(file, cv::IMREAD_GRAYSCALE);
	if(image.empty())
		return ;

	cv::imshow("gradient", image);
	cv::waitKey(0);

	wide = image.cols;
	height = image.rows;
	size = (wide+2) * (height+2);

	pixel = new int [size];
	float *Vx = new float [size];
	float *Vy = new float [size];
	createInputTexture();

	//create gradient
	for(int y = 0; y < height; y++)
		for(int x = 0; x < wide; x++)
		{
			float val1, val2;
			val1 = image.at<uchar>(height-1 - y, x);
			if(x+1 < wide)
			{ 
				val2 = image.at<uchar>(height-1 - y, x+1);
				Vx[IX(x+1, y+1)] = (val2 - val1) / 2;
			}
			else
				Vx[IX(x+1, y+1)] = 0;

			if(y+1 < height)
			{ 
				val2 = image.at<uchar>(height-2 - y, x);
				Vy[IX(x+1, y+1)] = (val2 - val1) / 2;
			}
			else
				Vy[IX(x+1, y+1)] = 0;
		}
	
	//test
	//getOutputTextureLIC(wide/2, height/2, Vx, Vy);
	//return;

	cv::Mat output = cv::Mat(height, wide, CV_8UC1);
	for(int y = 1; y <= height; y++)
		for(int x = 1; x <= wide; x++)
		{
			float val;
			if(DDA)
				val = getOutputTextureDDA(x, y, Vx, Vy);
			else
				val = getOutputTextureLIC(x, y, Vx, Vy);
			output.at<uchar>(height-y, x-1) = uchar(val * 255);
		}

	cv::imshow("Output", output);
	cv::waitKey(0);

	delete [] Vx;
	delete [] Vy;
}

