#include "NoQwt.h"

#include <QPen>
#include <QBrush>
#include <QString>
#include <QVector>
#include <QPainter>
#include <QLineF>
#include <QStaticText>
#include <QGraphicsScene>
#include <QWheelEvent>
#include <QStack>
#include <QDebug>

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtGraphicsView::Impl
{
    Impl() : scaleFactor(1.15) {}
    QStack<QPointF> scales;
    qreal scaleFactor;
};

NoQwtGraphicsView::NoQwtGraphicsView(QWidget* parent)
    : QGraphicsView(parent), pimpl(new Impl)
{
    setScene(new QGraphicsScene);

    setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
    setTransform(QTransform::fromScale(1,-1));
    //connect(view, SIGNAL(resizeEvent()), SLOG(scale_view()));
}

void NoQwtGraphicsView::wheelEvent(QWheelEvent* event)
{
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

    if(event->delta() > 0) { // Zoom in
        if(pimpl->scales.size() < 1000) {
            scale(pimpl->scaleFactor, pimpl->scaleFactor);
            pimpl->scales.push(QPointF(pimpl->scaleFactor, pimpl->scaleFactor));
        }
    } else { // Zooming out
        if(not pimpl->scales.isEmpty()) {
            const QPointF s = pimpl->scales.pop();
            scale(1.0 / s.x(), 1.0 / s.y());
        }
    }
}

void NoQwtGraphicsView::mouseReleaseEvent(QMouseEvent* event)
{
    if(event->button()==Qt::RightButton)
    {
        while(not pimpl->scales.isEmpty()) {
            const QPointF s = pimpl->scales.pop();
            scale(1.0 / s.x(), 1.0 / s.y());
        }
        event->accept();
    }
    else
        event->ignore();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtPlotCurve::Impl
{
    QPen pen;
    QBrush brush;
    QString title;
    QVector<QPointF> points;
    QRectF boundingRect;

    Impl(const QPen& p, const QBrush& b, const QString& t) : pen(p), brush(b), title(t)
    {
        //pen.setCosmetic(true);
    }
};

NoQwtPlotCurve::NoQwtPlotCurve(NoQwtPlot* parent, const QString& t, const QPen& p, const QBrush& b, const QString& tag)
    : QGraphicsObject(parent),
      pimpl(new Impl(p,b,t))
{
    parent->add_curve(this, tag);
}

NoQwtPlotCurve::~NoQwtPlotCurve()
{
    delete pimpl;
}

void NoQwtPlotCurve::addData(double x, double y)
{
    if(pimpl->points.isEmpty())
    {
        pimpl->boundingRect = QRectF(x, y, 0., 0.);
        pimpl->points.append(QPointF(x,y));
    }
    else if(not pimpl->boundingRect.contains(x, y))
    {
        prepareGeometryChange();
        pimpl->boundingRect.setLeft( qMin(pimpl->boundingRect.left(), x) );
        pimpl->boundingRect.setRight( qMax(pimpl->boundingRect.right(), x) );
        pimpl->boundingRect.setTop( qMin(pimpl->boundingRect.top(), y) );
        pimpl->boundingRect.setBottom( qMax(pimpl->boundingRect.bottom(), y) );
        if(pimpl->points.size() > 1)
            pimpl->points.append(pimpl->points.back());
        pimpl->points.append(QPointF(x,y));
        emit geometryChanged();
    }
    //update();
}

void NoQwtPlotCurve::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    painter->setPen(pimpl->pen);
    //painter->setBrush(pimpl->brush);
    painter->drawLines(pimpl->points);
}

QRectF NoQwtPlotCurve::boundingRect() const
{
    return pimpl->boundingRect;
}

void NoQwtPlotCurve::reset()
{
    prepareGeometryChange();
    pimpl->boundingRect = QRectF();
    pimpl->points.clear();
    emit geometryChanged();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtPlot::Impl
{
    QString title, x_title, y_title;
    QMap<QString,NoQwtPlotCurve*> curves;
    QRectF boundingRect;
    QPen axles_pen;

    Impl(const QString& t, const QString& x, const QString& y)
        : title(t), x_title(x), y_title(y)
    {
        axles_pen.setColor(Qt::black);
    }

    void paintGrid(QPainter *painter);
    void paintAxles(QPainter *painter);
    void paintAxleNocks(QPainter *painter);
};

void NoQwtPlot::Impl::paintGrid(QPainter*)
{
}

void NoQwtPlot::Impl::paintAxles(QPainter* painter)
{
    painter->setPen(axles_pen);
    painter->drawLine(boundingRect.topLeft(), boundingRect.topRight());
    painter->drawLine(boundingRect.topLeft(), boundingRect.bottomLeft());
}

void NoQwtPlot::Impl::paintAxleNocks(QPainter* painter)
{
    painter->setPen(axles_pen);

    const qreal box_w = boundingRect.right() - boundingRect.left();
    const qreal box_h = boundingRect.right() - boundingRect.left();

    qreal nock_step = 1.;
    while(box_w < 6 * nock_step)
        nock_step /= 10;

    qreal nock_pos = nock_step * int(boundingRect.left() / nock_step);
    do {
        const QPointF m_pos(nock_pos, boundingRect.top() + box_h/20);
        painter->drawLine(nock_pos, boundingRect.top(), m_pos.x(), m_pos.y());
    //    painter->drawText(m_pos, QString::number(nock_pos));
        nock_pos += nock_step;
    } while(nock_pos < boundingRect.right());
}

//-----------------------------------------------------------------------

NoQwtPlot::NoQwtPlot(QGraphicsScene* scene, const QString& t, const QString& x, const QString& y)
    : pimpl(new Impl(t,x,y))
{
    scene->addItem(this);
}

NoQwtPlot::~NoQwtPlot()
{
    delete pimpl;
}

QRectF NoQwtPlot::boundingRect() const
{
    return pimpl->boundingRect;
}

void NoQwtPlot::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    pimpl->paintGrid(painter);
    pimpl->paintAxles(painter);
    pimpl->paintAxleNocks(painter);
}

void NoQwtPlot::add_curve(NoQwtPlotCurve* curve, const QString& tag)
{
    pimpl->curves[tag] = curve;
    curve->setParentItem(this);
    connect(curve, SIGNAL(geometryChanged()), SLOT(childGeometryChanged()));
}

NoQwtPlotCurve* NoQwtPlot::curve(const QString& tag)
{
    return pimpl->curves.value(tag,NULL);
}

void NoQwtPlot::childGeometryChanged()
{
    prepareGeometryChange();

    pimpl->boundingRect = QRectF();
    foreach(const QGraphicsItem* c, childItems())
        pimpl->boundingRect |= c->boundingRect();
}

void NoQwtPlot::reset()
{
    foreach(NoQwtPlotCurve* c, pimpl->curves)
        c->reset();
}

