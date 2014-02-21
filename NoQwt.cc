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
    //setDragMode(ScrollHandDrag);
    //connect(view, SIGNAL(resizeEvent()), SLOG(scale_view()));
}

NoQwtGraphicsView::~NoQwtGraphicsView()
{
    delete pimpl;
}

void NoQwtGraphicsView::wheelEvent(QWheelEvent* event)
{
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

    if(event->delta() > 0) { // Zoom in
        if(pimpl->scales.size() < 1000) {
            scale(pimpl->scaleFactor, pimpl->scaleFactor);
            pimpl->scales.push(QPointF(pimpl->scaleFactor, pimpl->scaleFactor));
            emit scaled();
        }
    } else { // Zooming out
        if(not pimpl->scales.isEmpty()) {
            const QPointF s = pimpl->scales.pop();
            scale(1.0 / s.x(), 1.0 / s.y());
            emit scaled();
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
        emit scaled();
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
    Q_UNUSED(option);
    Q_UNUSED(widget);

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
    QVector<QLineF> v_grid, h_grid;
    QRectF boundingRect;
    QPen axles_pen, grid_pen;

    Impl(const QString& t, const QString& x, const QString& y)
        : title(t), x_title(x), y_title(y)
    {
        axles_pen.setColor(Qt::black);
        grid_pen.setColor(QColor(127,127,127,127));
        grid_pen.setStyle(Qt::DotLine);
    }
    QRectF view_box(const QGraphicsItem* i);

    void paintGrid(QPainter *painter);
    void paintAxles(QPainter *painter);
    void paintAxleNocks(QPainter *painter, const NoQwtPlot*);
};

void NoQwtPlot::Impl::paintGrid(QPainter* painter)
{
    painter->setPen(grid_pen);
    painter->drawLines(v_grid);
    painter->drawLines(h_grid);
}

void NoQwtPlot::Impl::paintAxles(QPainter* painter)
{
    painter->setPen(axles_pen);
    painter->drawLine(boundingRect.topLeft(), boundingRect.topRight());
    painter->drawLine(boundingRect.topLeft(), boundingRect.bottomLeft());
}

void NoQwtPlot::Impl::paintAxleNocks(QPainter* painter, const NoQwtPlot* plot)
{
    const QRectF box = view_box(plot);

    QPointF p;

    const QLineF left_box_line(box.topLeft(), box.bottomLeft());

    foreach(const QLineF& h, h_grid)
        if(h.intersect(left_box_line, &p))
            ;//painter->drawText(p, QString::number(h.y1()));

    const QLineF top_box_line(box.topLeft(), box.topRight());

    foreach(const QLineF& v, v_grid)
        if(v.intersect(top_box_line, &p))
            ;//painter->drawText(p, QString::number(v.x1()));
}

QRectF NoQwtPlot::Impl::view_box(const QGraphicsItem* i)
{
    const QGraphicsView* view = i->scene()->views().at(0);
    const QRect portRect = view->viewport()->rect();
    const QRectF sceneRect = view->mapToScene(portRect).boundingRect();
    return i->mapRectFromScene(sceneRect);
}

//-----------------------------------------------------------------------

NoQwtPlot::NoQwtPlot(QGraphicsScene* scene, const QString& t, const QString& x, const QString& y)
    : pimpl(new Impl(t,x,y))
{
    scene->addItem(this);
    QGraphicsView* view = scene->views().at(0);
    connect(view, SIGNAL(scaled()), this, SLOT(scaled()));

    setFlag(QGraphicsItem::ItemClipsChildrenToShape);
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
    Q_UNUSED(option);
    Q_UNUSED(widget);

    pimpl->paintGrid(painter);
    pimpl->paintAxles(painter);
    pimpl->paintAxleNocks(painter, this);
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
    scaled(); // не масштабировано, конечно, но boundingRect-то изменился.
}

void NoQwtPlot::reset()
{
    foreach(NoQwtPlotCurve* c, pimpl->curves)
        c->reset();
}

void NoQwtPlot::scaled()
{
    const QRectF box = pimpl->view_box(this);

    qreal nock_h_step = 10.;
    while(box.width() < 2 * nock_h_step)
        nock_h_step /= 10;

    pimpl->v_grid.clear();

    for(qreal nock_pos = 0.; nock_pos < pimpl->boundingRect.right(); nock_pos += nock_h_step)
        pimpl->v_grid.append(QLineF(nock_pos, pimpl->boundingRect.top(), nock_pos, pimpl->boundingRect.bottom()));
    for(qreal nock_pos = 0.; nock_pos > pimpl->boundingRect.left(); nock_pos -= nock_h_step)
        pimpl->v_grid.append(QLineF(nock_pos, pimpl->boundingRect.top(), nock_pos, pimpl->boundingRect.bottom()));

    qreal nock_v_step = 10.;
    while(box.height() < 2 * nock_v_step)
        nock_v_step /= 10;

    pimpl->h_grid.clear();

    for(qreal nock_pos = 0.; nock_pos < pimpl->boundingRect.bottom(); nock_pos += nock_v_step)
        pimpl->h_grid.append(QLineF(pimpl->boundingRect.left(), nock_pos, pimpl->boundingRect.right(), nock_pos));
    for(qreal nock_pos = 0.; nock_pos > pimpl->boundingRect.top(); nock_pos -= nock_v_step)
        pimpl->h_grid.append(QLineF(pimpl->boundingRect.left(), nock_pos, pimpl->boundingRect.right(), nock_pos));
}

