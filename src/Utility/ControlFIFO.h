#ifndef CONTROL_FIFO_H_
#define CONTROL_FIFO_H_

#include <map>
#include <string>
#include <tinyxml.h>


/*! \brief Callback for FIFO option handling */
class ControlCallback {
  public:
    //! \brief Callback on receiving a control block
    //! \param[in] context The XML block to parse
    virtual void OnControl(const TiXmlElement* context) = 0;

    //! \brief Get context name for callback
    virtual std::string GetContext() = 0;
};

/*! \brief This class enables simple app control over a FIFO.
 *  \details A fifo is opened, and users can write instructions to the fifo
 *           in XML format. These are then processed between time steps */
class ControlFIFO {
  public:
    //! \brief Default constructor
    ControlFIFO();

    //! \brief Default destructor. Tears down the opened fifo and removes the
    //         filesystem entry
    ~ControlFIFO();

    //! \brief Register a callback handler
    //! \param[in] callback The callback handler to register
    void registerCallback(ControlCallback& callback);

    //! \brief Open the fifo and prepare for receiving
    //! \param[in] name The name of the filesystem entry for the fifo
    bool open(const std::string& name="ifem-control");

    //! \brief Poll for new data in the fifo
    void poll();

  protected:
    std::string fifo_name; //!< The filesystem entry of our fifo
    int fifo;              //!< fifo handle
    typedef std::map<std::string, ControlCallback*> CallbackMap; //!< A map of callbacks
    CallbackMap callbacks; //!< Map of registered callbacks
};

#endif
