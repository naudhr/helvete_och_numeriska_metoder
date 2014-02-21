#ifndef __NO_QWT_H_
#define __NO_QWT_H_

#include <QGraphicsView>
#include <QGraphicsItem>
#include <QGraphicsObject>

//class QPen;
//class QBrush;
class QString;
class QPainter;
class QRectF;
class NoQwtPlotCurve;
class NoQwtPlotLegend;

class NoQwtGraphicsView : public QGraphicsView
{
    Q_OBJECT

    struct Impl;
    Impl* pimpl;

  public:
    ~NoQwtGraphicsView();
    NoQwtGraphicsView(QWidget* parent);
    void wheelEvent(QWheelEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);

  signals:
    void scaled();
};

class NoQwtPlot : public QGraphicsObject
{
    Q_OBJECT

    struct Impl;
    Impl* pimpl;

  public:
    NoQwtPlot(QGraphicsScene *parent, const QString& t, const QString& x, const QString& y);
    virtual ~NoQwtPlot();

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    virtual QRectF boundingRect() const;

    void addLegend();
    void removeLegend();
    NoQwtPlotLegend *legend();

    NoQwtPlotCurve* curve(const QString& label);
    void add_curve(NoQwtPlotCurve* curve, const QString& label);

  public slots:
    void reset();
    void scaled();

  private slots:
    void childGeometryChanged();
};

class NoQwtPlotCurve : public QGraphicsObject
{
    Q_OBJECT

    struct Impl;
    Impl* pimpl;

  public:

    NoQwtPlotCurve(NoQwtPlot *parent, const QString& t, const QPen& p, const QBrush& b, const QString& label);
    virtual ~NoQwtPlotCurve();

    const QPen& pen();
    const QBrush& brush();
    const QString& title();

    void addData(double x, double y);

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    virtual QRectF boundingRect() const;

  public slots:
    void reset();

  signals:
    void geometryChanged();
};
/*
class NoQwtPlotLegend : public QGraphicsObject
{
    Q_OBJECT

    struct Impl;
    Impl* pimpl;

public:
    explicit NoQwtPlotLegend(NoQwtPlot* parent);
    ~NoQwtPlotLegend();

    QPainterPath shape() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

protected slots:
    void itemRemove();
    void dataItemChange();
};
*/
#endif // GRAPHICSPLOTLEGEND_H
