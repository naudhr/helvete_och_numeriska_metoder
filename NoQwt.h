#ifndef __NO_QWT_H_
#define __NO_QWT_H_

#include <QGraphicsView>
#include <QGraphicsItem>
#include <QGraphicsObject>

//class QPen;
//class QBrush;
class QString;
//class QPainter;
class QRectF;
class NoQwtPlotCurve;
class NoQwtPlotLegend;
class QListWidgetItem;

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
    NoQwtPlotLegend* legend();

  signals:
    void scaled();
};

class NoQwtPlot : public QGraphicsObject
{
    Q_OBJECT

    struct Impl;
    Impl* pimpl;

  public:
    NoQwtPlot(NoQwtGraphicsView* parent, const QString& t, const QString& x, const QString& y);
    virtual ~NoQwtPlot();

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    virtual QRectF boundingRect() const;

    NoQwtPlotCurve* curve(const QString& label);
    void add_curve(NoQwtPlotCurve* curve, const QString& label);

  public slots:
    void reset();
    void scaled();
    void toggleVisibility(const QString& label);

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

    const QPen& pen() const;
    const QString& label() const;

    void addData(double x, double y);

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    virtual QRectF boundingRect() const;

  public slots:
    void reset();
    void toggleVisibility();

  signals:
    void geometryChanged();
};

class NoQwtPlotLegend : public QWidget
{
    Q_OBJECT

    struct Impl;
    Impl* pimpl;

  private slots:
    void someItemClicked(QListWidgetItem* item);

  protected:
    virtual QSize sizeHint() const;

  public:

    explicit NoQwtPlotLegend(NoQwtGraphicsView * parent);
    virtual ~NoQwtPlotLegend();

    void add_curve(const NoQwtPlotCurve* curve);
    void setVisibleSection(const QString& group, bool );

  signals:
    void toggleVisibility(const QString& label);
};

#endif // __NO_QWT_H_
