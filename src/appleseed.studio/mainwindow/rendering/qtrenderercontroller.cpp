
//
// This source file is part of appleseed.
// Visit https://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014-2018 Francois Beaune, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// Interface header.
#include "qtrenderercontroller.h"

// appleseed.renderer headers.
#include "renderer/api/log.h"

// Qt headers.
#include <QEventLoop>
#include <QIODevice>
#include <QMetaType>
#include <QString>
#include <QStringList>
#include <QTcpSocket>

using namespace renderer;

namespace appleseed {
namespace studio {

//
// QtRendererController class implementation.
//

Q_DECLARE_METATYPE(QAbstractSocket::SocketError);

namespace
{
    struct RegisterMetaTypes
    {
        RegisterMetaTypes()
        {
            qRegisterMetaType<QAbstractSocket::SocketError>();
        }
    };

    RegisterMetaTypes register_meta_types;
}

QtRendererController::QtRendererController()
  : m_event_loop(nullptr)
  , m_tcp_socket(nullptr)
{
    set_intention(ContinueRendering);
}

void QtRendererController::on_rendering_begin()
{
    set_intention(ContinueRendering);

    connect_to_notification_server();

    emit signal_rendering_begin();
}

void QtRendererController::on_rendering_success()
{
    disconnect_from_notification_server();

    emit signal_rendering_success();
}

void QtRendererController::on_rendering_abort()
{
    disconnect_from_notification_server();

    emit signal_rendering_abort();
}

void QtRendererController::on_rendering_pause()
{
    emit signal_rendering_pause();
}

void QtRendererController::on_rendering_resume()
{
    emit signal_rendering_resume();
}

void QtRendererController::on_frame_begin()
{
    const Intention intention = get_intention();

    if (intention == RestartRendering || intention == ReinitializeRendering)
        set_intention(ContinueRendering);

    emit signal_frame_begin();
}

void QtRendererController::on_frame_end()
{
    emit signal_frame_end();
}

void QtRendererController::on_progress()
{
    if (m_event_loop != nullptr)
        m_event_loop->processEvents();
}

void QtRendererController::set_intention(const Intention intention)
{
    m_intention = intention;
}

IRendererController::Intention QtRendererController::get_intention() const
{
    return m_intention;
}

void QtRendererController::slot_tcp_socket_connected()
{
    RENDERER_LOG_INFO("successfully connected to notification server.");
}

void QtRendererController::slot_tcp_socket_error(const QAbstractSocket::SocketError error)
{
    RENDERER_LOG_ERROR(
        "connection to notification server failed: %s",
        m_tcp_socket->errorString().toStdString().c_str());
}

void QtRendererController::slot_tcp_socket_ready_read()
{
    const QByteArray bytes = m_tcp_socket->readAll();
    const QString str = QString::fromAscii(bytes.data());
    const QStringList commands = str.split("\r\n", QString::SkipEmptyParts);

    RENDERER_LOG_DEBUG("received %d command%s.",
        commands.size(), commands.size() > 1 ? "s" : "");

    // todo: aggregate commands (e.g. reinitialize + restart -> reinitialize).
    for (const QString& command : commands)
    {
        if (command == "continue")
            set_intention(ContinueRendering);
        else if (command == "pause")
            set_intention(PauseRendering);
        else if (command == "terminate")
            set_intention(TerminateRendering);
        else if (command == "abort")
            set_intention(AbortRendering);
        else if (command == "restart")
            set_intention(RestartRendering);
        else if (command == "reinitialize")
            set_intention(ReinitializeRendering);
        else
        {
            // Ignore unknown commands.
            RENDERER_LOG_WARNING("received unknown command: \"%s\"", command.toStdString().c_str());
        }
    }
}

void QtRendererController::connect_to_notification_server()
{
    assert(
        (m_event_loop == nullptr && m_tcp_socket == nullptr) ||
        (m_event_loop != nullptr && m_tcp_socket != nullptr));

    if (m_event_loop != nullptr)
        return;

    RENDERER_LOG_INFO("connecting to notification server...");

    // Cannot parent `m_event_loop` to `this` since they typically belong to different threads.
    m_event_loop = new QEventLoop();

    // Cannot parent `m_tcp_socket` to `this` since they typically belong to different threads.
    m_tcp_socket = new QTcpSocket();

    connect(
        m_tcp_socket, SIGNAL(connected()),
        this, SLOT(slot_tcp_socket_connected()));

    connect(
        m_tcp_socket, SIGNAL(error(QAbstractSocket::SocketError)),
        this, SLOT(slot_tcp_socket_error(QAbstractSocket::SocketError)));

    connect(
        m_tcp_socket, SIGNAL(readyRead()),
        this, SLOT(slot_tcp_socket_ready_read()));

    m_tcp_socket->connectToHost("localhost", 54322, QIODevice::ReadOnly);
}

void QtRendererController::disconnect_from_notification_server()
{
    RENDERER_LOG_INFO("disconnecting from notification server...");

    if (m_tcp_socket != nullptr)
        m_tcp_socket->close();

    delete m_tcp_socket;
    m_tcp_socket = nullptr;

    delete m_event_loop;
    m_event_loop = nullptr;

    assert(m_event_loop == nullptr);
    assert(m_tcp_socket == nullptr);

    RENDERER_LOG_INFO("disconnected from notification server.");
}

}   // namespace studio
}   // namespace appleseed
