#include "NoQwt.h"

#include <QPen>
#include <QBrush>
#include <QString>
#include <QVector>
#include <QPainter>
#include <QLineF>
#include <QStaticText>
#include <QGraphicsScene>
#include <QGraphicsTextItem>
#include <QWheelEvent>
#include <QStack>
#include <QHBoxLayout>
#include <QListWidget>
#include <QPixmap>
#include <QLabel>
#include <QMessageBox>
#include <QDebug>

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtGraphicsView::Impl
{
    Impl() : scaleFactor(1.15) {}
    QStack<QPointF> scales;
    qreal scaleFactor;
    NoQwtPlotLegend* legend;
    //QVector<QPointF> h_nocks, v_nocks;
};

NoQwtGraphicsView::NoQwtGraphicsView(QWidget* parent)
    : QGraphicsView(parent), pimpl(new Impl)
{
    setScene(new QGraphicsScene);

    setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
    setTransform(QTransform::fromScale(1,-1));
    setDragMode(ScrollHandDrag);
    //connect(view, SIGNAL(resizeEvent()), SLOG(scale_view()));

    pimpl->legend = new NoQwtPlotLegend(this);
}

NoQwtGraphicsView::~NoQwtGraphicsView()
{
    delete pimpl;
}

NoQwtPlotLegend* NoQwtGraphicsView::legend()
{
    return pimpl->legend;
}

void NoQwtGraphicsView::wheelEvent(QWheelEvent* event)
{
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);

    if(event->delta() > 0) { // Zoom in
        if(pimpl->scales.size() < 100) {
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
        fitInView(scene()->items(Qt::AscendingOrder).at(0));
        pimpl->scales.clear();
        emit scaled();
    }
    else
        QGraphicsView::mouseReleaseEvent(event);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtPlotCurve::Impl
{
    QPen pen;
    QString title;
    QVector<QPointF> points;
    QRectF boundingRect;

    Impl(const QPen& p, const QString& t) : pen(p), title(t)
    {
        //pen.setCosmetic(true);
    }
};

const QPen& NoQwtPlotCurve::pen() const
{
    return pimpl->pen;
}

const QString& NoQwtPlotCurve::label() const
{
    return pimpl->title;
}

NoQwtPlotCurve::NoQwtPlotCurve(NoQwtPlot* parent, const QString& t, const QPen& p, const QString& tag)
    : QGraphicsObject(parent),
      pimpl(new Impl(p,t))
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
        pimpl->points.append(QPointF(x,y));
        emit geometryChanged();
    }
}

void NoQwtPlotCurve::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    Q_UNUSED(option);
    Q_UNUSED(widget);

    painter->setPen(pimpl->pen);
    painter->drawPolyline(pimpl->points.data(), pimpl->points.size());
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

void NoQwtPlotCurve::toggleVisibility()
{
    setVisible(not isVisible());
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtPlot::Impl
{
    NoQwtGraphicsView* view;
    QString title, x_title, y_title;
    QMap<QString,NoQwtPlotCurve*> curves;
    QVector<QLineF> v_grid, h_grid;
    QVector<QGraphicsTextItem*> labels;
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
    void paintAxleNocks(QPainter *painter);
};

void NoQwtPlot::Impl::paintGrid(QPainter* painter)
{
    painter->setPen(grid_pen);
    painter->drawLines(v_grid);
    painter->drawLines(h_grid);
}

void NoQwtPlot::Impl::paintAxles(QPainter* painter)
{
    //painter->setPen(axles_pen);
    //painter->drawLine(boundingRect.topLeft(), boundingRect.topRight());
    //painter->drawLine(boundingRect.topLeft(), boundingRect.bottomLeft());
}

void NoQwtPlot::Impl::paintAxleNocks(QPainter* painter)
{
    Q_UNUSED(painter);
}

QRectF NoQwtPlot::Impl::view_box(const QGraphicsItem* i)
{
    const QRect portRect = view->viewport()->rect();
    const QRectF sceneRect = view->mapToScene(portRect).boundingRect();
    return i->mapRectFromScene(sceneRect);
}

//-----------------------------------------------------------------------

NoQwtPlot::NoQwtPlot(NoQwtGraphicsView* view, const QString& t, const QString& x, const QString& y)
    : pimpl(new Impl(t,x,y))
{
    pimpl->view = view;
    view->scene()->addItem(this);
    connect(view, SIGNAL(scaled()), this, SLOT(scaled()));
    connect(view->legend(), SIGNAL(toggleVisibility(QString)), SLOT(toggleVisibility(QString)));

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
    pimpl->paintAxleNocks(painter);
}

void NoQwtPlot::add_curve(NoQwtPlotCurve* curve, const QString& tag)
{
    pimpl->curves[tag] = curve;
    curve->setParentItem(this);
    connect(curve, SIGNAL(geometryChanged()), SLOT(childGeometryChanged()));
    pimpl->view->legend()->add_curve(curve);
}

NoQwtPlotCurve* NoQwtPlot::curve(const QString& tag)
{
    return pimpl->curves.value(tag,NULL);
}

void NoQwtPlot::childGeometryChanged()
{
    prepareGeometryChange();
    pimpl->boundingRect = QRectF();
    foreach(const QGraphicsItem* c, pimpl->curves)
        pimpl->boundingRect |= c->boundingRect();
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
    while(box.width() < 4 * nock_h_step)
    {
        nock_h_step /= 2;
        if(box.width() < 2 * nock_h_step)
            nock_h_step /= 2;
        if(box.width() < 2 * nock_h_step)
            nock_h_step /= 2.5;
    }

    pimpl->v_grid.clear();

    for(qreal nock_pos = 0.; nock_pos < pimpl->boundingRect.right(); nock_pos += nock_h_step)
        pimpl->v_grid.append(QLineF(nock_pos, pimpl->boundingRect.top(), nock_pos, pimpl->boundingRect.bottom()));
    for(qreal nock_pos = 0.; nock_pos > pimpl->boundingRect.left(); nock_pos -= nock_h_step)
        pimpl->v_grid.append(QLineF(nock_pos, pimpl->boundingRect.top(), nock_pos, pimpl->boundingRect.bottom()));

    qreal nock_v_step = 10.;
    while(box.height() < 2 * nock_v_step)
    {
        nock_v_step /= 2;
        if(box.height() < 2 * nock_v_step)
            nock_v_step /= 2;
        if(box.height() < 2 * nock_v_step)
            nock_v_step /= 2.5;
    }

    pimpl->h_grid.clear();

    for(qreal nock_pos = 0.; nock_pos < pimpl->boundingRect.bottom(); nock_pos += nock_v_step)
        pimpl->h_grid.append(QLineF(pimpl->boundingRect.left(), nock_pos, pimpl->boundingRect.right(), nock_pos));
    for(qreal nock_pos = 0.; nock_pos > pimpl->boundingRect.top(); nock_pos -= nock_v_step)
        pimpl->h_grid.append(QLineF(pimpl->boundingRect.left(), nock_pos, pimpl->boundingRect.right(), nock_pos));

    foreach(QGraphicsTextItem* t, pimpl->labels)
        delete t;
    pimpl->labels.clear();

    QPointF p;
    foreach(const QLineF& h, pimpl->h_grid)
        foreach(const QLineF& v, pimpl->v_grid)
            if(h.intersect(v, &p))
            {
                QGraphicsTextItem* t = new QGraphicsTextItem(QString("%1 x %2").arg(p.x()).arg(p.y()), this);
                t->setFlag(ItemIgnoresTransformations);
                t->setPos(p);
                pimpl->labels.append(t);
            }
}

void NoQwtPlot::toggleVisibility(const QString& label)
{
    foreach(NoQwtPlotCurve* c, pimpl->curves)
        if(c->label() == label)
            c->toggleVisibility();
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

struct NoQwtPlotLegend::Impl
{
    QMap<QString, QListWidget*> groups;
};

//-----------------------------------------------------------------------

NoQwtPlotLegend::NoQwtPlotLegend(NoQwtGraphicsView* parent) : QWidget(parent), pimpl(new Impl)
{
    new QHBoxLayout(this);
    //QHBoxLayout* l = new QHBoxLayout(this);
    //l->addWidget(new QLabel("Legend:",this),1);
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
}

NoQwtPlotLegend::~NoQwtPlotLegend()
{
    delete pimpl;
}

QSize NoQwtPlotLegend::sizeHint() const
{
    return QSize(320, 150);
}

void NoQwtPlotLegend::add_curve(const NoQwtPlotCurve* curve)
{
    const QString group = curve->label().section(' ', -1);
    QListWidget* list = pimpl->groups.value(group,NULL);
    if(not list)
    {
        QHBoxLayout* grid = static_cast<QHBoxLayout*>(layout());
        grid->addWidget(list = new QListWidget(this));
        list->setSelectionMode(QAbstractItemView::NoSelection);
        list->addItem(group);
        pimpl->groups.insert(group,list);
        connect(list, SIGNAL(itemClicked(QListWidgetItem*)), SLOT(someItemClicked(QListWidgetItem*)));
    }
    QPixmap p(10,10);
    p.fill(curve->pen().color());
    new QListWidgetItem(p, curve->label().section(' ', -2, -2), list);
}

void NoQwtPlotLegend::setVisibleSection(const QString& group, bool v)
{
    if(QListWidget* list = pimpl->groups.value(group,NULL))
        list->setVisible(v);
    setVisible(true);
    move(parentWidget()->width() - width(), 0);
}

void NoQwtPlotLegend::someItemClicked(QListWidgetItem* item)
{
    foreach(QListWidget* list, pimpl->groups)
        if(list->row(item) >= 0)
        {
            emit toggleVisibility(item->text() + " " + list->item(0)->text());
            QFont f = item->font();
            f.setItalic(not f.italic());
            item->setFont(f);
        }
}
