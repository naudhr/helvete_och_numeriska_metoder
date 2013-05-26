#ifndef __CALCULUS_H_
#define __CALCULUS_H_

#include <QObject>
#include <QVector>
#include "params.h"
//struct Params;
//struct AnswerItem;

class CalculusEiler : public QObject
{
    Q_OBJECT
  signals:
    void a_step_done();
  public:
    QVector<AnswerItem> doWork(const Params& p);
};

class CalculusTrapeze : public QObject
{
    Q_OBJECT
  signals:
    void a_step_done();
  public:
    QVector<AnswerItem> doWork(const Params& p);
};

#endif // __CALCULUS_H_
