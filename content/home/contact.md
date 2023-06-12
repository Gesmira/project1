---
# An instance of the Contact widget.
widget: contact

# This file represents a page section.
headless: true

# Order that this section appears on the page.
weight: 130

title: Contact
subtitle:

content:
  email: gesmiramolla@gmail.com
  email: gmolla@nygenome.org
  address:
    city: New York
    region: NY
  # Automatically link email and phone
  autolink: true

  # Email form provider
  form:
    provider: netlify
    formspree:
      id:
    netlify:
      # Enable CAPTCHA challenge to reduce spam?
      captcha: false

  # Contact details (edit or remove options as required)


design:
  columns: '2'
---
