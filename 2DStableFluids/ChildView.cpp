
// ChildView.cpp : implementation of the CChildView class
//

#include "stdafx.h"
#include "2DStableFluids.h"
#include "ChildView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CChildView

CChildView::CChildView()
{
	windowSize = 600;
	dx = windowSize/fluidSolver.n;
	showDensity = true;
	showVelocity = false;
	showGrid = false;
	m_timer = 0;
	leftButton = false;
	rightButton = false;
	blackHot = true;

	fillDensity = 50.0f;
	maxFillDensity = 5000.0f;
	minFillDensity = 50.0f;
	deltaFillDensity = 50.0f;
}

CChildView::~CChildView()
{
}


BEGIN_MESSAGE_MAP(CChildView, CWnd)
	ON_WM_PAINT()
	ON_WM_TIMER()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_KEYDOWN()
END_MESSAGE_MAP()



// CChildView message handlers

BOOL CChildView::PreCreateWindow(CREATESTRUCT& cs) 
{
	if (!CWnd::PreCreateWindow(cs))
		return FALSE;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), NULL);

	return TRUE;
}

void CChildView::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	
	CRect rect;
	
	this->GetClientRect(rect);
    // double buffering to remove flickering
	
    CDC MemDC; 
    CBitmap MemBitmap;
    MemDC.CreateCompatibleDC(NULL);
    MemBitmap.CreateCompatibleBitmap(&dc,windowSize+1,windowSize+1);
    CBitmap *pOldBit=MemDC.SelectObject(&MemBitmap);


	COLORREF backgroundColor = blackHot ? RGB(0, 0, 0) : RGB(255, 255, 255);
	MemDC.FillSolidRect(0, 0, windowSize + 1, windowSize + 1, backgroundColor);

    // clear the screen
    CBrush brush;
    // repaint the window with solid black background
    brush.CreateSolidBrush(RGB(255,255,255));

	COLORREF qCircleColor = RGB(0,0,0);
    CPen qCirclePen(PS_SOLID, 5, qCircleColor);
    CPen* pqOrigPen = MemDC.SelectObject(&qCirclePen);

	int grid_number = fluidSolver.n;
	
	if (showDensity)
	{
		for (int cell_i = 0; cell_i < grid_number; cell_i++)
			for (int cell_j = 0; cell_j < grid_number; cell_j++)
			{
				int index = cell_i + cell_j * grid_number;
				double densityVal = fluidSolver.density[index];
				int intensity = (int)((1 - densityVal) * 255);

				if (intensity < 0) intensity = 0;
				if (intensity > 255) intensity = 255;

				if (blackHot)
				{
					// White circles on black background
					qCircleColor = RGB(intensity, intensity, intensity);
				}
				else
				{
					// Black circles on white background
					intensity = 255 - intensity;
					qCircleColor = RGB(intensity, intensity, intensity);
				}

				MemDC.FillSolidRect(cell_i * dx, cell_j * dx, dx, dx, qCircleColor);
			}
	}
	// Density color map
	else
	{
		// min and max pressure to normalize
		double minD = 0;
		double maxD = 1.5;
		

		for (int cell_i = 0; cell_i < grid_number; cell_i++) {
			for (int cell_j = 0; cell_j < grid_number; cell_j++) {
				int index = cell_i + cell_j * grid_number;
				double d = fluidSolver.density[index];
				double normD = (d - minD) / (maxD - minD);

				// Interpolate between colors using the normalized pressure
				COLORREF color = GetInterpolatedPressureColor(normD);
				MemDC.FillSolidRect(cell_i * dx, cell_j * dx, dx, dx, color);
			}
		}
	}

	if (showVelocity)
	{
		COLORREF qLineColor;
		if (showDensity)
		{
			// default to blue
			qLineColor = RGB(0, 0, 255);
		}
		else
		{
			// White if the density heat map is on for visibility
			qLineColor = RGB(255, 255, 255);
		}
		//draw new velocity
		CPen qLinePen(PS_SOLID, 1, qLineColor);
		MemDC.SelectObject(&qLinePen);
		vec2 *v;

		for (int cell_i = 0; cell_i < grid_number; cell_i++)
			for (int cell_j = 0; cell_j < grid_number; cell_j++)
			{
				v = fluidSolver.v(cell_i, cell_j);
				MemDC.MoveTo((cell_i) * dx, (cell_j) * dx);
				MemDC.LineTo((int) ((cell_i) * dx + v->x * dx) , (int) ( (cell_j) * dx + v->y * dx));
			}
	}

	//Draw Grid lines
	if (showGrid)
	{
		// Primal edges
		COLORREF qLineColor = RGB(125,125,125);
		CPen qLinePen(PS_SOLID, 1, qLineColor);
		MemDC.SelectObject(&qLinePen);
		for (int cell = 0; cell < grid_number; cell++)
		{
			// hrizental grid lines
			MemDC.MoveTo(0, (cell + 1) * dx);
			MemDC.LineTo(windowSize, (cell + 1) * dx);
			//vertical grid lines
			MemDC.MoveTo((cell + 1) * dx, 0);
			MemDC.LineTo((cell + 1) * dx, windowSize);
		}
	}

	dc.BitBlt(0,0,windowSize+1,windowSize+1,&MemDC,0,0,SRCCOPY);
    MemBitmap.DeleteObject();
    MemDC.DeleteDC();


	int TextWidth = 250;
	int TextHeight = 600;
	CDC MemDC1; 
    CBitmap MemBitmap1;
    MemDC1.CreateCompatibleDC(NULL);
    MemBitmap1.CreateCompatibleBitmap(&dc,TextWidth,TextHeight);
    CBitmap *pOldBit1=MemDC1.SelectObject(&MemBitmap1);
	//Background color for text
	MemDC1.FillSolidRect(windowSize+1, 0, TextWidth, TextHeight, RGB(0,20,20));
  
	MemDC1.SetTextColor(RGB(0,255,255));
	CString s1,s2;
	s1 = "n = ";
	s2.Format(_T("%d"),grid_number);
	s2 = s1+s2;
	MemDC1.TextOutW(3,10,s2);

	CString s3, s4;
	s3 = "Viscosity = ";
	s4.Format(_T("%.2f"), fluidSolver.viscosity);
	s3 = s3 + s4;
	MemDC1.TextOutW(3, 30, s3);

	CString s5, s6;
	s5 = "Buoyancy = ";
	s6.Format(_T("%.2f"), fluidSolver.buoyancy);
	s6 = s5 + s6;
	MemDC1.TextOutW(3, 50, s6);

	CString s7, s8;
	s7 = "Fill Density = ";
	s8.Format(_T("%.2f"), fillDensity);
	s8 = s7 + s8;
	MemDC1.TextOutW(3, 70, s8);

	MemDC1.SetTextColor(RGB(255,255,255));
	int row = 100;
	MemDC1.TextOutW(3,row,_T("User Interface Guide:"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("Z : Start Animation"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("Left button: inject smoke"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("Right button: stir fluid"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("R : Reset"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("V : Show/Hide Velocity"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("D : Show/Hide Density Heatmap"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("G : Show/Hide Gridlines"));
	row += 20;
	MemDC1.TextOutW(8, row, _T("B: Switch White Hot/Black Hot"));
	row += 20;
	MemDC1.TextOutW(8,row,_T("A : One more time step"));
	row += 20;
	MemDC1.TextOutW(8, row, _T("-/+: Decrease/Increase Viscosity"));
	row += 20;
	MemDC1.TextOutW(8, row, _T("Down/Up Arrow: Decr/Incr Buoyancy"));
	row += 20;
	MemDC1.TextOutW(8, row, _T("O/P Arrow: Decr/Incr Fill Tool Density"));


	dc.BitBlt(windowSize+1,0,TextWidth,TextHeight,&MemDC1,0,0,NOTSRCCOPY);
    MemBitmap1.DeleteObject();
    MemDC1.DeleteDC();

}

COLORREF CChildView::GetInterpolatedPressureColor(double value)
{
	// Clamp between 0 and 1
	if (value < 0.0) value = 0.0;
	if (value > 1.0) value = 1.0;

	// Define pressure bands and colors
	struct ColorStop {
		double position;
		COLORREF color;
	};

	static ColorStop stops[] = {
		{ 0.0, RGB(0, 0, 255) },     // Blue (low)
		{ 0.25, RGB(0, 255, 0) },    // Green (low-med)
		{ 0.5, RGB(255, 255, 0) },   // Yellow (med)
		{ 0.75, RGB(255, 165, 0) },  // Orange (med-high)
		{ 1.0, RGB(255, 0, 0) }      // Red (high)
	};

	// Find the two surrounding color stops
	for (int i = 0; i < 4; ++i) {
		if (value >= stops[i].position && value <= stops[i + 1].position) {
			double range = stops[i + 1].position - stops[i].position;
			double t = (value - stops[i].position) / range;

			BYTE r1 = GetRValue(stops[i].color);
			BYTE g1 = GetGValue(stops[i].color);
			BYTE b1 = GetBValue(stops[i].color);

			BYTE r2 = GetRValue(stops[i + 1].color);
			BYTE g2 = GetGValue(stops[i + 1].color);
			BYTE b2 = GetBValue(stops[i + 1].color);

			BYTE r = (BYTE)(r1 + t * (r2 - r1));
			BYTE g = (BYTE)(g1 + t * (g2 - g1));
			BYTE b = (BYTE)(b1 + t * (b2 - b1));

			return RGB(r, g, b);
		}
	}

	return RGB(255, 255, 255); // fallback white
}

void CChildView::OnTimer(UINT_PTR nIDEvent)
{
	if (leftButton) {
		int index = Find_Cell_Index(current_point);	 
		// Inject density
		fluidSolver.density_source[index] = fillDensity * fluidSolver.h;
	}

	fluidSolver.update();

	Invalidate(false);
	CWnd::OnTimer(nIDEvent);
}


void CChildView::OnLButtonDown(UINT nFlags, CPoint point)
{
	leftButton = true;
	current_point = old_point = point;
	CWnd::OnLButtonDown(nFlags, point);
}


void CChildView::OnLButtonUp(UINT nFlags, CPoint point)
{
	leftButton = false;
	CWnd::OnLButtonUp(nFlags, point);
}


void CChildView::OnRButtonDown(UINT nFlags, CPoint point)
{
	rightButton = true;
	current_point = old_point = point;
	CWnd::OnRButtonDown(nFlags, point);
}


void CChildView::OnRButtonUp(UINT nFlags, CPoint point)
{
	rightButton = false;

	CWnd::OnRButtonUp(nFlags, point);
}


void CChildView::OnMouseMove(UINT nFlags, CPoint point)
{
	if (leftButton || rightButton) {
		old_point = current_point;
		current_point = point;
	}

	if (rightButton) {
		int index = Find_Cell_Index(old_point);
		//Modify velocity
		fluidSolver.velocity_source[index] = vec2(current_point.x-old_point.x,current_point.y-old_point.y)*50.;
	}

	CWnd::OnMouseMove(nFlags, point);
}


void CChildView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{

	switch( nChar)
	{
	case 'a':
	case 'A':
		fluidSolver.update();
		break;
	case 'Z':
	case 'z':
		if (m_timer) {
			KillTimer(m_timer);
			m_timer = 0;
		}
		else
		{
		   m_timer = SetTimer(1, (int) (25), NULL);
		}
		break;
	case 'r':
	case 'R':
		fluidSolver.reset();
		Invalidate(false);
		break;
	case 'D':
	case 'd':
		showDensity = !showDensity;
		Invalidate(false);
		break;
	case 'v':
	case 'V':
		showVelocity = !showVelocity;
		Invalidate(false);
		break;
	case 'G': // Show/Hide the grid
	case 'g':
		showGrid = !showGrid;
		InvalidateRect(NULL,FALSE);
		break;
	case VK_OEM_MINUS:
		fluidSolver.viscosity -= fluidSolver.delta_viscosity;

		if (fluidSolver.viscosity < 0.0f)
		{
			fluidSolver.viscosity = 0.0f;
		}
		break;

	case VK_OEM_PLUS:
		fluidSolver.viscosity += fluidSolver.delta_viscosity;

		if (fluidSolver.viscosity > 10.0f)
		{
			fluidSolver.viscosity = 10.0f;
		}
		break;
	case VK_DOWN:
		fluidSolver.buoyancy -= fluidSolver.delta_buoyancy;

		if (fluidSolver.buoyancy < -5.0f)
		{
			fluidSolver.buoyancy = -5.0f;
		}
		break;

	case VK_UP:
		fluidSolver.buoyancy += fluidSolver.delta_buoyancy;

		if (fluidSolver.buoyancy > 5.0f)
		{
			fluidSolver.buoyancy = 5.0f;
		}
		break;

	case 'B':
	case 'b':
		blackHot = !blackHot;
		break;

	case 'o':
	case 'O':
		fillDensity -= deltaFillDensity;

		if (fillDensity < minFillDensity)
		{
			fillDensity = minFillDensity;
		}
		break;

	case 'p':
	case 'P':
		fillDensity += deltaFillDensity;

		if (fillDensity > maxFillDensity)
		{
			fillDensity = maxFillDensity;
		}
		break;
	}

	

	CWnd::OnKeyDown(nChar, nRepCnt, nFlags);
}

int CChildView::Find_Cell_Index(CPoint point)
{
	int x = point.x;
	int y = point.y;
	// set boundaries for mouse input to be inside the rectangular
	if (x >= windowSize)
		x = windowSize - 1;
	if (x < 0)
		x = 0;
	if (y >= windowSize)
		y = windowSize - 1;
	if (y < 0)
		y = 0;
	// finwindowSize the cell inwindowSizeex of mouse input
	int cell_i = x / dx;
	int cell_j = y / dx;

	return cell_i + fluidSolver.n * cell_j;
}
