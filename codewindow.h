#ifndef CODEWINDOW_H
#define CODEWINDOW_H

#include <QTextEdit>

class CodeWindow : public QTextEdit
{
    Q_OBJECT
public:
    explicit CodeWindow(QObject *parent = 0);

signals:

public slots:

};

#endif // CODEWINDOW_H
