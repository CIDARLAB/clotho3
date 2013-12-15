/* librecv: defines a handler for each channel
 * Needs to match up with
 *     `org.clothocad.core.aspects.Communicator.SendChannels`
 */

var librecv = new Object();
librecv._dispatch = new Object();
var commandMethods = new Object();

/* Concrete implementation of libtransport.receive() */
libtransport.receive = function (channel, message, flags) {
    if (channel in librecv._dispatch) {
        librecv._dispatch[channel](message);
    } else {
        alert("librecv.receive got bad message: channel " + channel + " is not defined");
    }
};


librecv._dispatch.showCommandResult = function (message) {
	alert("show command called?")
};

librecv._dispatch.showQueryCompletions = function (message) {
    var json = JSON.parse(message);
    command_bar.showSuggestions(json.suggestions);
};

librecv._dispatch.showTabList = function(message) {
    var tabs_list = JSON.parse(message);
    navigation_bar.showTabs(tabs_list);
};

librecv._dispatch.login = function(message) {
    libcookie.set("auth_key", message);
};

librecv._dispatch.logout = function(message) {
    libcookie.del("auth_key");
};

/**
 *When the dispatch is a commandList, the msgObj is comething
 *like:
 *
 
 {
    collectArgs
    sayArgs
    collectArgs
    updateArgs
 }
 
 Where the individual command's arguments are wrapped in an
object that will make sense onto to that function.  The particular
type of command in question is held as the 'command' field of the object

 */
librecv._dispatch.commandList = function(message) {
	//alert("command array received called?")
    var msgArray = JSON.parse(message);
    for (var i = 0; i < msgArray.length; i++) {
	    var cmdObj = msgArray[i];
        var cmd = cmdObj.command;
		//alert("command array: "+cmd)
        if (cmd in commandMethods) {
            commandMethods[cmd](cmdObj);
        } else {
            alert("commandList bad message: channel " + cmdObj.command);
        }
        
    }
}


commandMethods.callback = function (callbackMsg) {
	var command = "clotho.callback('"+callbackMsg.dooID+"');";
	var message = JSON.stringify({"command":command, "query":command});
	libsend.call("submitCommand", message);
};

commandMethods.inject = function (injectMsg) {
};

commandMethods.collect = function (collectMsg) {
    try {
        var data = collectMsg.data;
        for (var i = 0; i < data.length; i++) {
            collector.add(data[i]);
        }
//        alert("collect called, collector has x elements: " + Object.keys(data).length);
    } catch(err) {
        alert(err);
    }
};



commandMethods.say = function (message) {
    command_bar.postArticle(
        command_bar.wrap_dialog_bubble(
            command_bar.wrap_dialog_block('Server Response'),
            'dialog_response',
            'leftArrow'
        )
    );
};

commandMethods.alert = function (message) {
    alert(message);
};



commandMethods.showWidget = function (showWidgeMsg) {
    widget_space.show(showWidgeMsg);
};

commandMethods.update = function (updateMsg) {
    /* TODO - check if need to append javascript... may not need here but may in showView
    $("head").append("<script type='text/javascript'>"+js_string+"</script>");
    */
    widget_space.update(updateMsg);
};

commandMethods.removeWidget = function(remWidgeMsg) {
    //passes a view_id
    $("div[view_id="+remWidgeMsg+"]").remove();
}

commandMethods.moveWidget = function(moveWidgeMsg) {
    //passes a view_id
    $("div[view_id="+moveWidgeMsg+"]").remove();
}

commandMethods.addPage = function (addPageMsg) {
//    alert('addpage!' + addPageMsg.mode + '  ' + addPageMsg.ephemeral_link_page_id);
    libpage.add(addPageMsg.mode, addPageMsg.ephemeral_link_page_id);
};

commandMethods.removePage = function (remPageMsg) {
    window.close();
};


