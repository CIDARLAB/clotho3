'use strict';


Application.Extensions.service('SimpleService', function() {
    return {
        text : function() {return "Some Text From the Service"}
    }
});