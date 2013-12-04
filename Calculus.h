#ifndef __CALCULUS_H_
#define __CALCULUS_H_

#include <QObject>
#include <QVector>
#include "params.h"
//struct Params;
//struct AnswerItem;

class Calculus : public QObject
{
    Q_OBJECT
  signals:
    void a_step_done();
  public:
    void emit_a_step_done() {  emit a_step_done();  }
    virtual QVector<AnswerItem> doWork(const Params& p) = 0;
    virtual ~Calculus() {}
};

struct CalculusEiler : public Calculus
{
    QVector<AnswerItem> doWork(const Params& p);
};

struct CalculusTrapeze : public Calculus
{
    QVector<AnswerItem> doWork(const Params& p);
};

#endif // __CALCULUS_H_
