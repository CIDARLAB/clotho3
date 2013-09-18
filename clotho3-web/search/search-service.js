'use strict';

Application.Search.service('Searchbar', ['Clotho', 'ClientAPI', '$timeout', '$q', '$rootScope', function(Clotho, ClientAPI, $timeout, $q, $rootScope) {

    /******* config ******/
    var options = {};
    options.dateFilter = 'mediumTime';
    options.timeFilter = 'timestamp';

    /******* data *******/
    var log = {};

    var autocomplete = {};
    autocomplete.autocompletions = [];
    autocomplete.autoDetail = {};
    autocomplete.detailTemplate = {};
    autocomplete.detailModel = {};
    autocomplete.detailUUID = -1;

    //demo data
    log.entries = [
        {
            "text" : "Welcome to Clotho!",
            "from" : "server",
            "class" : "success",
            "timestamp" : Date.now()
        }
    ];


    /****** display ******/
    var display = {};
    display.query = '';
    display.queryHistory = [];
    display.autocomplete = false; // autocomplete list
    display.autocompleteDetail = false; //pane to left of autocomplete
    display.autocompleteDetailInfo = false; // e.g. command or author
    display.help = false; // help menu far right
    display.log = false; // activity log
    display.logSnippet = false; // snippet right of log button

    // todo - DRY
    display.genLogPos = function() {
        var target = document.getElementById('searchbar_logbutton');
        display.logpos = {
            left : (target.offsetLeft + (target.scrollWidth / 2) - 200) + "px",
            top : (target.offsetTop + target.scrollHeight)  + "px"
        };
    };

    display.genAutocompletePos = function() {
        var target = document.getElementById('searchBarInput');
        display.autocompletePos = {
            left : (target.offsetLeft) + "px",
            top : (target.offsetTop + target.scrollHeight + 6)  + "px"
        };
    };

    display.show = function (field) {
        if (!display[field])
            display[field] = true;
    };

    display.hide = function(field) {
        if (display[field])
            display[field] = false;
    };

    display.toggle = function(field) {
        display[field] = !display[field];
    };

    display.detail = function(uuid) {
        if (typeof uuid == 'undefined') return;

        if (uuid != autocomplete.detailUUID) {
            display.hide('autocompleteDetailInfo');
            autocomplete.detailUUID = uuid;
        }

        Clotho.autocompleteDetail(autocomplete.detailUUID).then(function(result) {
            autocomplete.autoDetail = result;
            display.show('autocompleteDetail');
        });
    };

    display.undetail = function() {
        display.hide('autocompleteDetail');
        display.hide('autocompleteDetailInfo');
        autocomplete.detailModel = {};
    };

    //todo - avoid using index in case sort - have to namespace
    display.detailInfo = function (type, index) {
        //choose template
        switch (type) {
            case 'command' : {
                autocomplete.detailTemplate = 'search/detail-command.html';
                break;
            }
            case 'author' : {
                autocomplete.detailTemplate = 'search/detail-author.html';
                break;
            }
            default : {}
        }
        //choose model
        autocomplete.detailModel = autocomplete.autoDetail.sharables[index];
        if (type == "author")
            autocomplete.detailModel = autocomplete.detailModel.author;

        display.show('autocompleteDetailInfo');
    };


    /***** functions *****/

    function receiveMessage (data) {
        log.unread = (!!log.unread && !display.log) ? log.unread + 1 : 1;
        log.entries.unshift(data);
        display.show('logSnippet');
        //todo - cancel if new request comes in
        $timeout( function() {
            display.hide('logSnippet');
        }, 5000);

    }

    var execute = function (command) {
        console.log("this would be run: " + command);
        display.hide('autocomplete');
        display.undetail();
    };

    var submit = function (query) {
        if (!query)
            query = display.query;

        console.log(query);


        if (!!query) {
            var submission = {class : 'info', from : 'client', text: query};
            display.queryHistory.push(submission);
            ClientAPI.say(submission);
            Clotho.submit(query).then(function(result){
                display.query = '';
                ClientAPI.say({text: result || 'Command received'});
            });
            //display.autocomplete = false;
            display.undetail();
        }
    };

    /****** listeners *****/

    Clotho.listen("activityLog", function (data) {
        receiveMessage(data);
    }, 'searchbar');
    
    return {
        options : options,
        display : display,
        log : log,
        setQuery : function(item, $event) {
            if (typeof $event != 'undefined')
                $event.preventDefault();
            /*if (item.type != 'command') {
                display.undetail();
            }*/
            display.query = !!item.value ? item.value : item.text;
        },
        autocomplete : autocomplete,
        submit : submit,
        execute : execute
    }

}]);
