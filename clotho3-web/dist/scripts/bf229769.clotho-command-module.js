angular.module("clotho.tokenizer", []), angular.module("clotho.commandbar", ["clotho.foundation", "clotho.tokenizer", "ui.keypress"]), angular.module("ui.keypress", []).factory("keypressHelper", ["$parse",
  function(a) {
    var b = {
      8: "backspace",
      9: "tab",
      13: "enter",
      27: "esc",
      32: "space",
      33: "pageup",
      34: "pagedown",
      35: "end",
      36: "home",
      37: "left",
      38: "up",
      39: "right",
      40: "down",
      45: "insert",
      46: "delete"
    }, c = function(a) {
      return a.charAt(0).toUpperCase() + a.slice(1)
    };
    return function(d, e, f, g) {
      var h, i = [];
      h = e.$eval(g["ui" + c(d)]), angular.forEach(h, function(b, c) {
        var d, e;
        e = a(b), angular.forEach(c.split(" "), function(a) {
          d = {
            expression: e,
            keys: {}
          }, angular.forEach(a.split("-"), function(a) {
            d.keys[a] = !0
          }), i.push(d)
        })
      }), f.bind(d, function(a) {
        var c = !(!a.metaKey || a.ctrlKey),
          f = !! a.altKey,
          g = !! a.ctrlKey,
          h = !! a.shiftKey,
          j = a.keyCode;
        "keypress" === d && !h && j >= 97 && 122 >= j && (j -= 32), angular.forEach(i, function(d) {
          var i = d.keys[b[j]] || d.keys[j.toString()],
            k = !! d.keys.meta,
            l = !! d.keys.alt,
            m = !! d.keys.ctrl,
            n = !! d.keys.shift;
          i && k === c && l === f && m === g && n === h && e.$apply(function() {
            d.expression(e, {
              $event: a
            })
          })
        })
      })
    }
  }
]), angular.module("ui.keypress").directive("uiKeydown", ["keypressHelper",
  function(a) {
    return {
      link: function(b, c, d) {
        a("keydown", b, c, d)
      }
    }
  }
]), angular.module("ui.keypress").directive("uiKeypress", ["keypressHelper",
  function(a) {
    return {
      link: function(b, c, d) {
        a("keypress", b, c, d)
      }
    }
  }
]), angular.module("ui.keypress").directive("uiKeyup", ["keypressHelper",
  function(a) {
    return {
      link: function(b, c, d) {
        a("keyup", b, c, d)
      }
    }
  }
]), angular.module("clotho.tokenizer").directive("autoGrow", ["$timeout",
  function(a) {
    return {
      require: "?ngModel",
      link: function(b, c, d, e) {
        var f = (c.css("paddingLeft"), c.css("paddingRight"), 100),
          g = angular.element("<span></span>").css({
            position: "absolute",
            top: "-10000px",
            left: "-10000px",
            fontSize: c.css("fontSize"),
            fontFamily: c.css("fontFamily"),
            "white-space": "pre"
          });
        c.after(g);
        var h = function() {
          var a = c.val().replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/&/g, "&amp;");
          g.html("" !== a ? a : c[0].placeholder);
          var b = g[0].offsetWidth + 26,
            d = Math.max(b, f) + "px";
          c.css("width", d)
        };
        e ? b.$watch(function() {
          return e.$viewValue
        }, h) : c.bind("keyup keydown blur", h), a(h)
      }
    }
  }
]), angular.module("clotho.commandbar").service("CommandBar", ["Clotho", "ClientAPI", "ClothoCommandHistory", "Debug", "ClothoSchemas", "$timeout", "$q", "$document",
  function(a, b, c, d, e, f, g, h) {
    var i = (new d("Command Bar", "#ffbb55"), function() {
        return angular.element(h[0].querySelector("clotho-command-bar"))
      }),
      j = function() {
        return angular.element(h[0].querySelector("clotho-command-bar [clotho-tokenizer]"))
      }, k = function() {
        return angular.element(h[0].querySelector("clotho-command-bar [clotho-reference-autocomplete]"))
      }, l = "query",
      m = function() {
        k().focus()
      }, n = function(a) {
        k().scope().setQueryString(a)
      }, o = {};
    return o.log = !1, o.logSnippet = !1, o.toggle = function(a, b) {
      o[a] = angular.isDefined(b) ? b : !o[a]
    }, o.toggleActivityLog = function() {
      o.log = !o.log, o.log && (log.unread = "")
    }, {
      display: o,
      setQuery: n,
      getCommandBarElement: i,
      getTokenizerElement: j,
      getCommandBarInput: k,
      commandBarInputModel: l,
      focusInput: m,
      setInput: n
    }
  }
]), angular.module("clotho.commandbar").directive("clothoCommandBar", ["Clotho", "CommandBar", "ClothoCommandHistory", "terminalAsideOptions", "$location", "$window", "$compile", "$timeout", "$clothoModal",
  function(a, b, c, d, e, f, g, h, i) {
    return {
      restrict: "EA",
      scope: !0,
      templateUrl: "views/_command/commandbar.html",
      controller: ["$scope", "$element", "$attrs",
        function(a) {
          a.logEntries = c.entries, a.display = b.display;
          var d = !1;
          a.toggleLogin = function(a) {
            d = angular.isDefined(a) ? a : !d, d ? i.create({
              title: "Clotho Login",
              "template-url": "'views/_command/simpleLogin.html'"
            }) : i.destroy()
          };
          var e = null,
            f = 1e4;
          a.startLogTimeout = function() {
            a.cancelLogTimeout(), e = h(function() {
              a.display.toggle("logSnippet", !1)
            }, f)
          }, a.cancelLogTimeout = function() {
            h.cancel(e)
          }, a.$watchCollection("logEntries", function() {
            a.display.toggle("logSnippet", !0), a.startLogTimeout()
          })
        }
      ],
      link: function(b) {
        b.toggleTerminalAside = c.toggleTerminal, b.showMeHow = function() {
          a.query({
            name: "Learning Clotho"
          }, {
            mute: !0
          }).then(function(b) {
            a.startTrail(b[0].id)
          })
        }, b.goHome = function() {
          e.path("/")
        }, b.aboutClotho = function() {
          e.path("/about")
        }, b.teamClotho = function() {
          e.path("/team")
        }
      }
    }
  }
]), angular.module("clotho.commandbar").controller("TerminalCtrl", ["$scope", "ClothoCommandHistory",
  function(a, b) {
    a.logEntries = b.entries
  }
]), angular.module("clotho.commandbar").value("terminalAsideOptions", {
  visible: !1
}).directive("terminalAside", ["$http", "$q", "$templateCache", "$window", "$animate", "$compile", "terminalAsideOptions", "ClothoCommandHistory", "Clotho",
  function(a, b, c, d, e, f, g, h, i) {
    return {
      restrict: "EA",
      templateUrl: "views/_command/terminalAside.html",
      replace: !0,
      scope: {
        title: "=?asideTitle",
        contentUrl: "=?asideContentUrl"
      },
      controller: ["$scope", "$element", "$attrs",
        function(a) {
          a.submit = function() {
            a.terminalQuery && i.submit(a.terminalQuery)
          }
        }
      ],
      link: function(d, e) {
        function j(d) {
          return b.when(c.get(d) || a.get(d)).then(function(a) {
            return angular.isObject(a) ? (c.put(d, a.data), a.data) : a
          })
        }
        d.$hide = function() {
          g.visible = !1
        }, d.$show = function() {
          g.visible = !0
        }, d.$toggle = function(a) {
          g.visible = angular.isDefined(a) ? a : !g.visible, angular.element("body").attr("aside-status", g.visible ? "active" : "")
        }, d.$watch("contentUrl", function(a) {
          a && j(a).then(function(a) {
            console.info("ASIDE TEMPLATE" + a), d.content = f(a)(d)
          })
        }), i.listen("toggleTerminalActive", d.$toggle, d), d.$watch(function() {
          return g.visible
        }, function(a) {
          a ? (e.addClass("active"), h.setLastView()) : e.removeClass("active")
        })
      }
    }
  }
]).directive("terminalAsideTrigger", ["Clotho", "$timeout", "terminalAsideOptions",
  function(a, b, c) {
    return {
      restrict: "A",
      template: '<div id="terminalAsideTrigger" ng-click="toggle()" ng-attr-status="{{activeClass ? \'active\' : \'\'}}" ng-class="{active : activeClass}"></div>',
      replace: !0,
      scope: !0,
      link: function(d) {
        d.toggle = function() {
          a.trigger("toggleTerminalActive"), b(function() {
            d.activeClass = c.visible, angular.element("body").attr("aside-status", c.visible ? "active" : "")
          })
        }
      }
    }
  }
]), angular.module("clotho.commandbar").directive("logEntries", function() {
  return {
    restrict: "A",
    templateUrl: "views/_command/logEntries.html",
    scope: {
      entries: "=logEntries"
    }
  }
}), angular.module("clotho.commandbar").controller("loginCtrl", ["$scope", "$timeout", "$location", "Clotho", "ClothoAuth", "PubSub", "Facebook",
  function(a, b, c, d, e, f, g) {
    function h() {
      a.cred.password = "", a.cred.confirm = ""
    }

    function i() {
      a.cred = {
        username: "",
        password: "",
        confirm: "",
        personId: ""
      }
    }
    a.createMode = !1, a.notification = {}, i(), f.on("auth:login", function() {
      c.replace(), c.path("/settings")
    }), a.login = function() {
      d.login(a.cred.username, a.cred.password).then(function(b) {
        b ? a.notification = {
          "class": "alert-success",
          message: "Log in Success"
        } : (a.notification = {
          "class": "alert-danger",
          message: "Log in Error"
        }, h())
      }, function() {})
    }, a.$watch("cred.username", function(b) {
      b && d.get(b, {
        mute: !0
      }).then(function(b) {
        b && b.id ? (a.retrieved = b, a.cred.personId = b.id, a.createMode = !1) : (a.retrieved = null, a.cred.personId = "")
      }, function() {
        a.retrieved = null, a.cred.personId = ""
      })
    }), a.logout = function() {
      g.logout().then(function() {
        i(), a.createMode = !1, a.notification = {}, a.showFacebookLogin = !0, a.retrieved = null, a.facebookRetrieved = !1
      })
    }, a.facebookLogin = function() {
      g.login().then(function(b) {
        a.cred.username = b.email, a.facebookRetrieved = !0, d.get(b.email, {
          mute: !0
        }).then(function(c) {
          if (angular.isEmpty(c)) {
            console.log("get was null, creating");
            var e = g.convertToPersonSharable(b);
            d.create(e).then(function(c) {
              a.notification = {
                "class": "alert-success",
                message: "Account added to Clotho! Choose a password"
              }, a.retrieved = e, a.cred.username = b.email, a.cred.personId = c
            }, function() {
              a.notification = {
                "class": "alert-success",
                message: "You've imported that account already! Choose a password"
              }, a.retrieved = e, a.cred.username = b.email, a.cred.personId = e.id
            }).then(function() {
              a.createMode = !0
            })
          } else console.log("get returned", c), a.notification = {
            "class": "alert-success",
            message: "Enter your password to login!"
          }, a.cred.username = b.email, a.cred.personId = b.email
        })
      }, function() {
        a.notification = {
          "class": "alert-danger",
          message: "Logging into facebook went wrong..."
        }
      })
    }, a.createAccount = function() {
      d.createUser(a.cred.username, a.cred.password).then(function(b) {
        console.log("create user?", b), a.notification = b ? {
          "class": "alert-success",
          message: "User " + a.cred.username + "created!"
        } : {
          "class": "alert-error",
          message: "Account creation unsuccessful"
        }
      }, function(b) {
        a.notification = {
          "class": "alert-error",
          message: "Error Creating... check console"
        }, console.error("account creation error", b)
      })
    }
  }
]), angular.module("clotho.commandbar").directive("clothoFunctionExecutor", ["$filter", "Clotho", "ClothoSchemas",
  function(a, b, c) {
    return {
      scope: {
        "function": "=",
        onExecute: "&?"
      },
      templateUrl: "views/_command/executor.html",
      link: function(d) {
        function e() {
          d.functionArgs = {}
        }

        function f(a, b) {
          var c = [];
          return angular.forEach(a.args, function(a) {
            c.push(b[a.name])
          }), c
        }
        d.isPrimitiveField = c.isPrimitiveField, d.schemaReadable = c.mapSchemaIdToName, d.capitalize = function(a) {
          return a.substring(0, 1).toUpperCase() + a.substr(1)
        }, d.queryWrapper = function(d, e) {
          return b.autocomplete(d).then(function(b) {
            return angular.isUndefined(e) ? b : a("filter")(b, function(a) {
              return "function" == e ? c.isFunction(a) : c.isInstanceOfSchema(a, e)
            })
          })
        }, d.executeFunction = function() {
          b.run(d.function.id, f(d.function, d.functionArgs)).then(function(a) {
            d.functionResult = a, d.onExecute({
              $result: a
            })
          })
        }, d.clearArguments = function() {
          e(), d.functionResult = null
        }, d.setArgument = function(a, c) {
          b.get(c, {
            mute: !0
          }).then(function(b) {
            d.functionArgs[a] = b
          })
        }, d.$watch("function.id", function() {
          e()
        })
      }
    }
  }
]), angular.module("clotho.tokenizer").directive("clothoAutocompleteListing", function() {
  return {
    restrict: "EA",
    scope: {
      autocompletions: "=",
      query: "=",
      active: "=",
      hasFocus: "=",
      triggerHide: "=",
      select: "&",
      passedPlacement: "@?",
      forceVisible: "@?"
    },
    replace: !0,
    templateUrl: "views/_command/autocompleteListing.html",
    link: function(a) {
      a.isVisible = function() {
        return a.forceVisible === !0 || a.forceVisible === !1 ? a.forceVisible : a.hasFocus && !a.triggerHide && a.autocompletions.length
      }, a.isActive = function(b) {
        return a.active == b
      }, a.selectActive = function(b) {
        a.active = b
      }, a.selectMatch = function(b) {
        a.select({
          activeIdx: b
        })
      }
    }
  }
}), angular.module("clotho.tokenizer").directive("clothoAutocompleteMatch", ["ClothoSchemas",
  function(a) {
    return {
      restrict: "EA",
      replace: !0,
      scope: {
        index: "=",
        match: "=",
        query: "=",
        passedPlacement: "@"
      },
      templateUrl: "views/_command/autocompleteMatch.html",
      link: function(b, c, d) {
        b.$watch(function() {
          return d.active
        }, function(a) {
          b.active = b.$eval(a)
        }), b.$watch("match", function(c) {
          b.iconClass = a.determineSharableIcon(a.dirtyDetermineType(c))
        })
      }
    }
  }
]).filter("clothoAutocompleteHighlight", function() {
  function a(a) {
    return a.replace(/([.?*+^$[\]\\(){}|-])/g, "\\$1")
  }
  return function(b, c) {
    return c ? b.replace(new RegExp(a(c), "gi"), "<strong>$&</strong>") : b
  }
}), angular.module("clotho.tokenizer").directive("clothoToken", ["$parse", "$document", "Clotho", "clothoTokenFactory", "ClothoSchemas",
  function(a, b, c, d, e) {
    return {
      restrict: "E",
      templateUrl: "views/_command/token.html",
      scope: {
        token: "=?",
        tokenId: "@?",
        tokenModel: "=?",
        tokenName: "@?",
        tokenActivePass: "@?",
        popupPlacement: "@?",
        popupTrigger: "@?",
        onClick: "&?",
        onRemove: "&?"
      },
      link: function(a, f) {
        a.$watch("token.model", function(b) {
          var c = e.dirtyDetermineType(b);
          a.labelClass = "label-" + e.typeToColorClass(c), a.iconClass = e.determineSharableIcon(c), angular.isDefined(a.token) && a.token.isSharable() && (a.labelClass += " isSharable")
        }), a.$watch("tokenModel", function(b) {
          angular.isEmpty(b) || (a.token = new d(b))
        }), a.$watch("tokenId", function(b) {
          angular.isEmpty(b) || (a.tokenName = b, c.get(b, {
            mute: !0
          }).then(function(b) {
            angular.isEmpty(b) || (a.token = new d(b))
          }))
        }), a.$watch("tokenActivePass", function(b) {
          angular.isDefined(b) && a.token.active(a.$eval(b))
        }), f.on("click", function(b) {
          a.onClick({
            $event: b,
            $token: a.token
          }), a.$apply(function() {
            a.token.active(!a.token.active())
          })
        });
        var g = function(b) {
          8 === b.which ? (b.preventDefault(), f.remove(), a.$destroy()) : 27 === b.which && (b.preventDefault(), a.$apply(function() {
            a.token.active(!1)
          }))
        };
        a.$watch("token.active()", function(a) {
          f.toggleClass("active", !! a), a ? b.on("keydown", g) : b.off("keydown", g)
        }), a.$on("$destroy", function() {
          a.onRemove({
            $token: a.token
          }), b.off("keydown", g)
        })
      }
    }
  }
]), angular.module("clotho.tokenizer").factory("clothoTokenCollectionFactory", ["clothoTokenFactory",
  function(a) {
    function b(a) {
      var b = this;
      this.tokens = [], (angular.isArray(a) || angular.isObject(a)) && angular.forEach(a, function(a) {
        b.addToken(a)
      }), c = function(a) {
        return angular.isDefined(a) ? b.tokens[a] : b.tokens[b.tokens.length - 1]
      }
    }
    var c = angular.noop;
    return b.prototype.addToken = function(b) {
      this.tokens.push(new a(b))
    }, b.prototype.hasLength = function() {
      return this.tokens.length > 0
    }, b.prototype.inRange = function(a) {
      return a > -1 && a < this.tokens.length
    }, b.prototype.getToken = function(a) {
      return c(a)
    }, b.prototype.indexOf = function(a) {
      return this.tokens.indexOf(a)
    }, b.prototype.removeToken = function(a) {
      return this.inRange(a) ? this.tokens.splice(a, 1) : !1
    }, b.prototype.removeAll = function() {
      this.tokens.length = 0
    }, b.prototype.removeActiveTokens = function() {
      for (var a = [], b = 0; b < this.tokens.length; b++)
        if (c(b).active()) {
          var d = this.removeToken(b);
          a.push(d)
        }
      return a
    }, b.prototype.setActive = function(a) {
      return this.inRange(a) ? (c(a).active(!0), !0) : !1
    }, b.prototype.unsetActive = function() {
      angular.forEach(this.tokens, function(a) {
        a.active(!1)
      })
    }, b.prototype.isActive = function(a) {
      if (angular.isDefined(a)) return c(a).active();
      for (var b = 0; b < this.tokens.length; b++)
        if (c(b).active()) return !0;
      return !1
    }, b.prototype.setLastActive = function() {
      return this.hasLength() && c().active(!0)
    }, b.prototype.isLastActive = function() {
      return this.hasLength() && c().active()
    }, b.prototype.clearLast = function() {
      return this.removeToken(this.tokens.length - 1)
    }, b
  }
]), angular.module("clotho.tokenizer").factory("clothoTokenFactory", ["$q", "Clotho",
  function(a, b) {
    function c(c) {
      var d = this;
      d.model = c, d.activeState = !1, d.fullSharablePromise = this.isSharable() ? b.get(d.model.id, {
        mute: !0
      }).then(function(a) {
        return d.fullSharable = a, a
      }) : a.when(null)
    }
    return c.prototype.readable = function() {
      return this.model.name || this.model
    }, c.prototype.id = function() {
      return this.model.id || null
    }, c.prototype.isAmbiguous = function() {
      return angular.isArray(this.model)
    }, c.prototype.isSharable = function() {
      return !this.isAmbiguous() && angular.isDefined(this.model.id)
    }, c.prototype.active = function(a) {
      return angular.isDefined(a) ? void(this.activeState = !! a) : this.activeState === !0
    }, c
  }
]), angular.module("clotho.tokenizer").directive("clothoReferenceListener", ["ClothoReferenceDelimiter", "$parse",
  function(a, b) {
    return function(c, d, e) {
      d.on("keydown", function(f) {
        event.which === a.keycode && (c.$eval(e.clothoReferencePreventDefault) && f.preventDefault(), b(e.clothoReferenceCallback)(c, {
          $event: f,
          $element: d,
          $delimiter: a.symbol
        }))
      })
    }
  }
]), angular.module("clotho.tokenizer").constant("ClothoReferenceDelimiter", {
  symbol: "@",
  keycode: 64
}).service("ClothoReference", ["ClothoReferenceDelimiter",
  function(a) {
    function b(b, d) {
      return angular.isUndefined(d) || !d.length ? "" : d.replace(c, function(c, d, e, f, g) {
        return d + b(e, a.symbol, f, g)
      })
    }
    var c = new RegExp("(^|[^a-zA-Z0-9-_\\.])" + a.symbol + "([A-Za-z0-9-_\\.@]+)", "gi"),
      d = "http://www.clothocad.org/data/";
    this.convert = b, this.convertHtml = angular.bind(null, b, function(a) {
      return '<clotho-token token-id="' + a + '" popup-trigger="mouseenter"></clotho-token>'
    }), this.convertMarkdown = angular.bind(null, b, function(a) {
      return "[" + a + "](" + d + a + ")"
    }), this.convertWiki = angular.bind(null, b, function(a) {
      return "[" + d + a + " " + a + "]"
    })
  }
]).directive("clothoReferenceParser", ["ClothoReference", "$compile", "$timeout", "$parse",
  function(a, b) {
    return function(c, d, e) {
      c.$watch(e.clothoReferenceParser, function(f) {
        "markdown" === e.clothoReferenceType ? d.text(a.convertMarkdown(f)) : "wiki" === e.clothoReferenceType ? d.text(a.convertWiki(f)) : (d.html(a.convertHtml(f)), b(d.contents())(c))
      })
    }
  }
]).filter("clothoReferenceParse", ["ClothoReference",
  function(a) {
    return function(b, c, d) {
      if (angular.isFunction(d)) return d(b);
      if (angular.isUndefined(c)) return b;
      switch (c) {
        case "html":
          return a.convertHtml(b);
        case "wiki":
          return a.convertWiki(b);
        case "markdown":
          return a.convertMarkdown(b)
      }
    }
  }
]), angular.module("clotho.commandbar").directive("clothoTerminalInput", ["Clotho", "ClientAPI",
  function(a, b) {
    return {
      restrict: "E",
      templateUrl: "views/_command/terminalInput.html",
      scope: {
        placeholder: "@",
        autocompleteTrigger: "@?"
      },
      link: function(c, d) {
        c.focusInput = function() {
          d.find("[clotho-reference-autocomplete]").focus()
        }, c.submit = function() {
          b.say({
            from: "client",
            channel: "submit",
            "class": "info",
            text: c.query
          }), a.submit(c.query).then(function(a) {
            b.say({
              from: "server",
              channel: "submit",
              "class": "success",
              text: a
            }), c.query = ""
          })
        }
      }
    }
  }
]), angular.module("clotho.tokenizer").directive("clothoReferenceAutocomplete", ["$q", "$parse", "$timeout", "$compile", "$filter", "$document", "Clotho", "ClothoReferenceDelimiter",
  function(a, b, c, d, e, f, g, h) {
    var i = [8, 9, 13, 27, 37, 38, 39, 40];
    return i.push(h.keycode), {
      restrict: "A",
      require: "?ngModel",
      scope: {
        query: "=ngModel",
        autocompleteTrigger: "=?",
        autocompleteTriggerInclude: "=?",
        autocompleteClearOnSelect: "=?",
        autocompleteClearOnEnter: "=?",
        forceVisible: "=?",
        autocompletions: "=?",
        autocompleteBlockInput: "=?",
        autocompleteFilter: "=?",
        autocompleteHasFocus: "=?",
        autocompleteOnSelect: "&?",
        autocompleteOnKeydown: "&?",
        autocompleteOnQuery: "&?",
        autocompleteOnBackout: "&?",
        autocompleteOnEnter: "&?",
        autocompleteDelimiter: "@?",
        autocompletePopupPosition: "@?",
        autocompleteWaitTime: "@?"
      },
      link: function(a, b) {
        function j(b) {
          if (a.autocompleteTrigger) {
            var c = new RegExp(".*" + h.symbol + "(.+?)(?=[\\s'\"" + (angular.isDefined(a.autocompleteDelimiter) ? a.autocompleteDelimiter : "") + "]|$)", "ig").exec(b);
            return angular.isEmpty(c) ? "" : c[1]
          }
          return b
        }

        function k() {
          a.query = ""
        }

        function l() {
          a.autocompletions = [], a.activeIdx = -1
        }

        function m() {
          l(), a.$apply(function() {
            a.autocompleteHasFocus = !1
          })
        }

        function n(a) {
          return s.test(a.charAt(0))
        }

        function o(a) {
          return {
            $query: a,
            $auto: j(a)
          }
        }

        function p() {
          a.autocompleteTrigger && (a.triggerHide = !0), a.query.length && a.autocompletions.length ? (l(), b[0].focus(), a.$digest()) : (b[0].blur(), m())
        }

        function q(d) {
          a.autocompleteHasFocus && a.activeIdx < 0 && (angular.isUndefined(d) || f[0].activeElement != b[0] && !b[0].contains(d.target)) && c(m)
        }
        var r = angular.copy(i);
        b.attr({
          autocomplete: "off",
          autocorrect: "off",
          autocapitalize: "off",
          spellcheck: "false"
        });
        var s = /^['"].*/,
          t = angular.element("<clotho-autocomplete-listing></clotho-autocomplete-listing>");
        t.attr({
          autocompletions: "autocompletions",
          active: "activeIdx",
          select: "select(activeIdx)",
          "has-focus": "autocompleteHasFocus",
          query: "query",
          "force-visible": "{{forceVisible}}",
          "trigger-hide": "triggerHide",
          "passed-placement": "{{autocompletePopupPosition}}"
        });
        var u = function(b) {
          s.test(b) && (b = b.substring(1)), 0 !== b.length && g.autocomplete(b).then(function(b) {
            b && b.length ? (a.activeIdx = -1, a.autocompleteFilter && (b = e("filter")(b, a.autocompleteFilter)), a.autocompletions = e("limitTo")(b, 10)) : l()
          })
        };
        a.setQueryString = function(b) {
          var c = angular.isEmpty(b) ? a.query : b;
          a.query = c, l()
        }, a.select = function(d) {
          var e = d > -1 ? a.autocompletions[d] : null,
            f = o(a.query);
          if (a.autocompleteOnSelect(angular.extend({
              $item: e
            }, f)), a.autocompleteClearOnSelect === !0) k();
          else {
            var g = "" + (a.autocompleteTriggerInclude === !0 ? "@" : "") + (e.id ? e.id : e);
            a.query = a.query.replace("@" + f.$auto, g)
          }
          c(function() {
            l(), b[0].focus()
          }, 0, !1)
        };
        var v;
        a.$watch("query", function(b, d) {
          var e = o(a.query),
            f = e.$auto;
          a.autocompleteOnQuery(angular.extend({
            $old: d
          }, e)), f.length ? (a.autocompleteHasFocus = !0, a.autocompleteWaitTime > 0 ? (v && c.cancel(v), v = c(function() {
            u(f)
          }, a.autocompleteWaitTime)) : u(f)) : l()
        }), a.$watch("autocompleteBlockInput", function(a) {
          b.toggleClass("input-disabled", !! a)
        }), a.autocompleteDelimiter && r.push(a.autocompleteDelimiter), a.$watch("autocompleteTrigger", function(b) {
          a.triggerHide = !! b
        }), b.bind("keydown", function(b) {
          if (a.autocompleteOnKeydown({
              $event: b,
              $keycode: b.which
            }), 50 === b.which && b.shiftKey === !0 && (a.autocompleteTrigger && (a.triggerHide = !1), a.$digest()), -1 === r.indexOf(b.which)) return void(a.autocompleteBlockInput && b.preventDefault());
          if (b.which === a.autocompleteDelimiter) "" == a.query || n(a.query) || a.$apply(function() {
            a.select(1 == a.autocompletions.length ? 0 : -1)
          });
          else if (8 === b.which) {
            if (a.query.length) return;
            a.$apply(function() {
              a.autocompleteOnBackout({
                $event: b
              })
            })
          } else if (40 === b.which) a.autocompletions.length && (a.activeIdx = (a.activeIdx + 1) % a.autocompletions.length, a.$digest());
          else if (38 === b.which) a.autocompletions.length && (a.activeIdx = (a.activeIdx ? a.activeIdx : a.autocompletions.length) - 1, a.$digest());
          else {
            if (37 === b.which) return;
            if (39 === b.which) return;
            13 === b.which || 9 === b.which ? (a.activeIdx >= 0 ? a.$apply(function() {
              a.select(a.activeIdx)
            }) : 13 == b.which && a.$apply(function() {
              a.autocompleteOnEnter(angular.extend({
                $event: b
              }, o(a.query))), a.autocompleteClearOnEnter === !0 && k()
            }), a.$digest()) : 27 === b.which && p(b)
          }
          b.preventDefault()
        });
        var w = function(a) {
          27 === a.which && p(a)
        };
        b.on("focus", function() {
          f.off("keydown", w), c(function() {
            a.autocompleteHasFocus = !0
          })
        }), b.on("blur", function() {
          a.autocompleteHasFocus && f.on("keydown", w)
        }), b.on("paste", function() {
          c(function() {
            a.setQueryString()
          })
        }), a.$on("$locationChangeSuccess", function() {
          setTimeout(m)
        }), f.bind("click", q), a.$on("$destroy", function() {
          f.unbind("click", q)
        }), a.autocompleteHasFocus = a.autocompleteHasFocus === !0, a.query = a.query || "", l(), b.after(d(t)(a))
      }
    }
  }
]), angular.module("clotho.tokenizer").directive("clothoReferenceTokenizer", ["$parse", "$timeout", "Clotho", "ClientAPI", "Debug", "clothoTokenCollectionFactory", "ClothoSchemas", "ClothoCommandHistory",
  function(a, b, c, d, e, f, g) {
    var h = new e("clothoReferenceTokenizer", "#ee9955");
    return {
      restrict: "E",
      templateUrl: "views/_command/clothoReferenceTokenizer.html",
      scope: {
        startingTags: "=?",
        passedPlaceholder: "@?placeholder"
      },
      controller: ["$scope", "$element", "$attrs",
        function(a, e) {
          function f() {
            a.query = ""
          }

          function i() {
            a.tokenCollection.removeAll()
          }
          this.addToken = function(b) {
            h.log("adding token", b), a.tokenCollection.addToken(b)
          }, this.removeToken = function(b) {
            h.log("removing token", b), a.tokenCollection.removeToken(b)
          }, this.handleBackout = function() {
            a.tokenCollection.isLastActive() ? a.tokenCollection.clearLast() : b(function() {
              a.tokenCollection.setLastActive()
            })
          }, this.onTokenActive = function() {}, this.handleSelect = function(a) {
            this.addToken(a)
          }, this.handleEnter = function(b, c) {
            b.preventDefault(), a.allowPrimitive === !0 && c.length ? this.addToken(c) : this.disambiguate(b)
          }, this.disambiguate = function() {
            var b = a.tokenCollection.tokens;
            if (b.length) {
              var e = [],
                h = b[0];
              if (angular.forEach(b, function(a) {
                  e.push(a.readable())
                }), d.say({
                  from: "client",
                  channel: "run",
                  "class": "info",
                  tokens: angular.copy(b),
                  text: "disambiguating: " + e.join(" ")
                }), g.isFunction(h.model)) {
                var j = angular.map(b.slice(1), function(a) {
                  return a.isSharable() ? a.model.id : a.model
                });
                c.run(h.model.id, j).then(function(a) {
                  d.say({
                    from: "server",
                    channel: "submit",
                    "class": "success",
                    text: a
                  })
                })
              } else g.isView(b[0]) ? console.log("need to handle view from command bar") : c.edit(b[0].id);
              f(), i()
            } else a.currentPlaceholder = "Must select from dropdown"
          }, this.focusInput = function() {
            e.find("[clotho-reference-autocomplete]").focus()
          }
        }
      ],
      controllerAs: "tokenCtrl",
      link: function(a) {
        function b() {
          a.currentFilter = null, a.currentPlaceholder = null
        }

        function c(b) {
          a.currentFilter = b
        }

        function d(b) {
          a.currentPlaceholder = b
        }

        function e(b) {
          a.allowPrimitive = angular.isDefined(b) ? !! b : !0
        }

        function h() {
          a.allowPrimitive = !1
        }

        function i() {
          c(function() {
            return !1
          })
        }

        function j() {
          i(), a.blockInput = !0
        }

        function k() {
          a.blockInput = !1, b()
        }
        a.tokenCollection = new f(a.startingTags), a.$watchCollection("tokenCollection.tokens", function(a) {
          if (a && a.length > 0) {
            var b = a[0],
              f = a.length;
            g.isFunction(b.model) ? b.fullSharablePromise.then(function(a) {
              if (a.args)
                if (a.args.length > f - 1) {
                  var b = a.args[f - 1],
                    l = g.isPrimitiveField(b.type);
                  l ? (i(), d("Enter " + b.type), e()) : (c(function(a) {
                    return g.isInstanceOfSchema(a, b.type)
                  }), d(b.type), h())
                } else j(), d("All arguments defined"), h();
              else k(), e()
            }) : (j(), d("Does not accept arguments"))
          } else k(), h()
        })
      }
    }
  }
]);
