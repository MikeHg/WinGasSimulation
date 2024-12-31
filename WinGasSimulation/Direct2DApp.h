#pragma once
#include <windows.h>
#include <d2d1.h>
#include <dwrite.h>

#pragma comment(lib, "d2d1.lib")
#pragma comment(lib, "dwrite.lib")

// A class to encapsulate all Direct2D related operations
// Draw box, particles etc.
class Direct2DApp {
public:
    Direct2DApp() : pFactory(nullptr), pRenderTarget(nullptr), pBrush(nullptr), pDWriteFactory(nullptr), pTextFormat(nullptr) {}

    ~Direct2DApp() {
        Cleanup();
    }

    // Initialize the necessary Direct2D and DirectWrite resources
    HRESULT Initialize(HWND hwnd) {
        HRESULT hr = D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, &pFactory);
        if (SUCCEEDED(hr)) {
            hr = DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(IDWriteFactory), reinterpret_cast<IUnknown**>(&pDWriteFactory));
        }
        if (SUCCEEDED(hr)) {
            hr = pDWriteFactory->CreateTextFormat(
                L"Arial",                // Font family
                nullptr,                 // Font collection
                DWRITE_FONT_WEIGHT_NORMAL,
                DWRITE_FONT_STYLE_NORMAL,
                DWRITE_FONT_STRETCH_NORMAL,
                18.0f,                   // Font size
                L"en-US",                // Locale
                &pTextFormat
            );
        }
        if (SUCCEEDED(hr)) {
            pTextFormat->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
            pTextFormat->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
        }
        if (SUCCEEDED(hr)) {
            // create a standard Maxwell distribution
            pFactory->CreatePathGeometry(&pPathGeometry);
            pPathGeometry->Open(&pSink);
            pSink->BeginFigure(
                D2D1::Point2F(box->getMaxwell(0, 0) + BOX_X * 2 + 40, box->getMaxwell(0, 1)),
                D2D1_FIGURE_BEGIN_HOLLOW
            );

            for (int i = 1; i < NUM_BOXES; i++) {
                pSink->AddBezier(
                    D2D1::BezierSegment(
                        D2D1::Point2F(box->getFirstCP(i - 1, 0) + BOX_X * 2 + 40, box->getFirstCP(i - 1, 1)),
                        D2D1::Point2F(box->getSecondCP(i - 1, 0) + BOX_X * 2 + 40, box->getSecondCP(i - 1, 1)),
                        D2D1::Point2F(box->getMaxwell(i, 0) + BOX_X * 2 + 40, box->getMaxwell(i, 1)))
                );
            }

            pSink->EndFigure(D2D1_FIGURE_END_OPEN);
            pSink->Close();
        }
        if (SUCCEEDED(hr)) {
            CreateGraphicsResources(hwnd);
        }
        return hr;
    }

    // Clean up the resources
    void Cleanup() {
        SafeRelease(&pPathGeometry);
        SafeRelease(&pBrush);
        SafeRelease(&pRenderTarget);
        SafeRelease(&pFactory);
        SafeRelease(&pDWriteFactory);
        SafeRelease(&pTextFormat);
    }

    // Handle the window paint event
    void OnPaint(HWND hwnd);

    // Handle window resizing
    void OnResize(HWND hwnd) {
        if (pRenderTarget != nullptr) {
            RECT rc;
            GetClientRect(hwnd, &rc);
            D2D1_SIZE_U size = D2D1::SizeU(rc.right, rc.bottom);
            pRenderTarget->Resize(size);
        }
    }

private:
    ID2D1Factory* pFactory;
    ID2D1HwndRenderTarget* pRenderTarget;
    ID2D1SolidColorBrush* pBrush;
    ID2D1SolidColorBrush* pBrush2;
    IDWriteFactory* pDWriteFactory;
    IDWriteTextFormat* pTextFormat;
    ID2D1GeometrySink* pSink = NULL;
    ID2D1PathGeometry* pPathGeometry;

    // Helper method to create necessary graphics resources
    void CreateGraphicsResources(HWND hwnd) {
        if (pRenderTarget == nullptr) {
            RECT rc;
            GetClientRect(hwnd, &rc);

            D2D1_SIZE_U size = D2D1::SizeU(rc.right, rc.bottom);

            pFactory->CreateHwndRenderTarget(
                D2D1::RenderTargetProperties(),
                D2D1::HwndRenderTargetProperties(hwnd, size),
                &pRenderTarget
            );

            pRenderTarget->CreateSolidColorBrush(
                D2D1::ColorF(D2D1::ColorF::Black),
                &pBrush
            );

            pRenderTarget->CreateSolidColorBrush(
                D2D1::ColorF(D2D1::ColorF::Red),
                &pBrush2
            );
        }
    }

    // Helper method to discard graphics resources when not needed
    void DiscardGraphicsResources() {
        SafeRelease(&pBrush);
        SafeRelease(&pRenderTarget);
    }

    // Template to safely release COM objects
    template <class T>
    void SafeRelease(T** ppT) {
        if (*ppT) {
            (*ppT)->Release();
            *ppT = nullptr;
        }
    }
};