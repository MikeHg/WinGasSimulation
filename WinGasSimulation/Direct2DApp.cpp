#include "pch.h"
#include "Direct2DApp.h"
#include <tchar.h>

#define BALLRADIUS 4

void Direct2DApp::OnPaint(HWND hwnd) {
    PAINTSTRUCT ps;
    TCHAR wszText[256];
    UINT32 cTextLength;
    D2D1_RECT_F layoutRect;
    
    BeginPaint(hwnd, &ps);

    CreateGraphicsResources(hwnd);

    pRenderTarget->BeginDraw();

    pRenderTarget->Clear(D2D1::ColorF(D2D1::ColorF::White));

    D2D1_RECT_F rectangle = D2D1::RectF(10 - BALLRADIUS, 10 - BALLRADIUS, 10+BOX_X*2 + BALLRADIUS, 10+BOX_Y*2 + BALLRADIUS);
    D2D1_ELLIPSE ellipse; 

    pRenderTarget->DrawRectangle(&rectangle, pBrush);

    ellipse = D2D1::Ellipse(D2D1::Point2F((float)(2.0L * box->getParticleLocation(0)[0] + 10.0), (float)(2.0L * box->getParticleLocation(0)[1] + 10.0))
        , BALLRADIUS+2, BALLRADIUS+2);
    pRenderTarget->FillEllipse(ellipse, pBrush2);

    for (int i = 1; i < NUM_PARTS; i++) {
        ellipse = D2D1::Ellipse(D2D1::Point2F((float)(2.0L*box->getParticleLocation(i)[0] + 10.0), (float)(2.0L *box->getParticleLocation(i)[1] + 10.0))
            , BALLRADIUS, BALLRADIUS);
        pRenderTarget->FillEllipse(ellipse, pBrush);
    }

 
    // all the text on the right hand side
    _stprintf_s(wszText, L"Step number: %u ", box->getStep());
    cTextLength = (UINT32)wcslen(wszText);
    layoutRect = D2D1::RectF(BOX_X * 2 + 40, 70, BOX_X * 2 + 550, 90);
    pRenderTarget->DrawText(
        wszText,
        cTextLength,
        pTextFormat,
        layoutRect,
        pBrush
    );


    if (box->getEnergyError() < 100.0L)
        _stprintf_s(wszText, L"Max. rel. error total energy: %f ", (double)box->getEnergyError());
    else
        _stprintf_s(wszText, L"Max. rel. error total energy: ---.------ ");

    cTextLength = (UINT32)wcslen(wszText);
    layoutRect = D2D1::RectF(BOX_X * 2 + 40, 10, BOX_X * 2 + 550, 30);

    pRenderTarget->DrawText(
        wszText,
        cTextLength,
        pTextFormat,
        layoutRect,
        pBrush
    );

    double vAbs = sqrt(box->getParticleVelocity(0)[0]* box->getParticleVelocity(0)[0] 
        + box->getParticleVelocity(0)[1] * box->getParticleVelocity(0)[1] 
        + box->getParticleVelocity(0)[2] * box->getParticleVelocity(0)[2]);
    _stprintf_s(wszText, L"Velocity (absolut) of particle 1: %f ", vAbs);
    cTextLength = (UINT32)wcslen(wszText);
    layoutRect = D2D1::RectF(BOX_X * 2 + 40, 40, BOX_X * 2 + 550, 60);

    pRenderTarget->DrawText(
        wszText,
        cTextLength,
        pTextFormat,
        layoutRect,
        pBrush
    );

    _stprintf_s(wszText, L"Velocity distribution of all particles");
    cTextLength = (UINT32)wcslen(wszText);
    layoutRect = D2D1::RectF(BOX_X * 2 + 40, 450, BOX_X * 2 + 550, 470);

    pRenderTarget->DrawText(
        wszText,
        cTextLength,
        pTextFormat,
        layoutRect,
        pBrush
    );

    // Draw the default Maxwell distribution
    pRenderTarget->DrawGeometry(pPathGeometry, pBrush, 1.f);

    for (int i = 0; i < NUM_BOXES; i++) {
        rectangle = D2D1::RectF(BOX_X * 2 + 40 + i*(BOX_X / NUM_BOXES), BOX_Y + 150.0, BOX_X * 2 + 40 + (i+1)*(BOX_X / NUM_BOXES), BOX_Y + 150.0 - box->getDistribV(i));
        pRenderTarget->DrawRectangle(&rectangle, pBrush);
    }

    HRESULT hr = pRenderTarget->EndDraw();
    if (hr == D2DERR_RECREATE_TARGET)
    {
        DiscardGraphicsResources();
    }

    EndPaint(hwnd, &ps);
}